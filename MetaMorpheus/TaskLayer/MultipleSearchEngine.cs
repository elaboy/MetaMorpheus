using Easy.Common.Extensions;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.Data;
using System.IO;
using System.Linq;
using System.Text.Json;
using System.Threading.Tasks;
using Proteomics.AminoAcidPolymer;
using static Nett.TomlObjectFactory;
using System.Drawing.Imaging;
using System.Reflection.Metadata.Ecma335;
using ICSharpCode.SharpZipLib;
using MathNet.Numerics;
using TorchSharp;

namespace TaskLayer
{
    public class MultipleSearchEngine
    {
        public List<List<Modification>> CombinationOfModifications { get; set; }
        public List<KeyValuePair<double, Modification[]>> CombinationsWithAddedMass { get; set; }
        public double[] MassArray { get; set; }
        public List<Tuple<double,Protein, MsDataScan, double, int>> ProteinListInferedFromGPTMD { get; private set; }
        /// <summary>
        /// double is MM score and int is precursorChargeState
        /// </summary>
        public List<Tuple<PeptideWithSetModifications, List<MatchedFragmentIon>, double, int>> peptideGroup { get; private set; }
        public MultipleSearchEngine(List<FilteredPsmTSV> psmList, List<Modification> listOfMods, int numberOfVariableMods, List<Modification> fixedMods, MsDataFile dataFile,
            bool allCombos)
        {
            List<List<Modification>> comboList = new();
            CombinationsWithAddedMass = new();
            Mc(listOfMods, ref comboList, numberOfVariableMods, allCombos);
            CombinationOfModifications = comboList;
            CombinationsWithAddedMass.Add(CombinationOfModifications.Select(x => new KeyValuePair<double, Modification[]>(
                key: x.Select(x => x.MonoisotopicMass.Value).Sum(), value: x.Select(x => x).ToArray())));
            CombinationsWithAddedMass = CombinationsWithAddedMass.OrderBy(x => x.Key).ToList();
            MassArray = CombinationsWithAddedMass.Select(x => x.Key).ToArray();
            SetProteinsInfered(psmList, dataFile);
            //SetPeptideGroup(fixedMods);
        }

        private void SetProteinsInfered(List<FilteredPsmTSV> psmList, MsDataFile dataFile)
        {
            ProteinListInferedFromGPTMD = new();
            Parallel.ForEach(psmList, psm =>
            {
                ProteinListInferedFromGPTMD.Add(new Tuple<double, Protein, MsDataScan, double, int>(double.Parse(psm.PrecursorMass),
                    new Protein(psm.BaseSeq, psm.ProteinAccession), dataFile.GetOneBasedScan(int.Parse(psm.ScanNumber)),
                    double.Parse(psm.Score), int.Parse(psm.Charge)));
            });
        }

        /// <summary>
        /// Reads the psmList and sets the peptide group to be used
        /// </summary>
        /// <param name="psmList"></param>
        //private void SetPeptideGroup(List<Modification> fixedMods)
        //{
        //    peptideGroup = new();
        //    foreach (var protein in ProteinListInferedFromGPTMD)
        //    {
        //        peptideGroup.Add(new Tuple<double, IEnumerable<PeptideWithSetModifications>>(protein.Item1,
        //            protein.Item2.Digest(new DigestionParams("top-down", 0), fixedMods, new List<Modification>())));
        //    }
        //}

        public void Run(string savingJsonPath, List<Modification> fixedMods, MsDataFile dataFile)
        {
            peptideGroup = new();

            foreach(var protein in ProteinListInferedFromGPTMD)
            {
                var fixedModdedPeptide = protein.Item2.Digest(new DigestionParams("top-down", 0),
                    fixedMods, new List<Modification>());


                var deltaMass = fixedModdedPeptide.First().MonoisotopicMass - protein.Item1;

                var possibleComboMods = GetCombinationsThatFitDelta(deltaMass);

                var products = new List<Product>();
                //Declare not to use list of mods position and mod

                foreach (var combo in possibleComboMods)
                {

                    var comboToTry = protein.Item2.Digest(new DigestionParams("top-down", 0),
                        fixedMods, combo);


                    //Exclude those that are in the no go list Check
                    foreach (var peptide in comboToTry)
                    {
                        peptide.Fragment(protein.Item3.DissociationType ?? DissociationType.HCD,
                            FragmentationTerminus.Both, products);

                        var match = MetaMorpheusEngine.MatchFragmentIons(new Ms2ScanWithSpecificMass(protein.Item3,
                                protein.Item3.SelectedIonMZ.Value, protein.Item5,
                                dataFile.FilePath, new CommonParameters()),
                            products, new CommonParameters());

                        if(match.Count > 0)
                            peptideGroup.Add(
                            new Tuple<PeptideWithSetModifications, List<MatchedFragmentIon>, double, int>( peptide, match, protein.Item4, protein.Item5));

                        //Add to the no go list generated from not improving position and mod
                    }
                }
            }


            SaveAsJson(savingJsonPath);
        }

        /// <summary>
        /// Returns the best candidate for the variant peptides.
        /// Not the one with the highest count of matched fragment ions but the one with the best coverage.
        /// </summary>
        /// <param name="peptideVariants"></param>
        /// <param name="modCombo"></param>
        /// <param name="possibleCombosThatFitDelta"></param>
        /// <param name="dataScan"></param>
        /// <param name="dataFilePath"></param>
        /// <param name="precursorChargeState"></param>
        /// <param name="products"></param>
        /// <returns></returns>
        private Tuple<PeptideWithSetModifications, List<MatchedFragmentIon>, double, int> GetBestVariant(
            IEnumerable<PeptideWithSetModifications> peptideVariants, List<Modification> modCombo,
            List<List<Modification>> possibleCombosThatFitDelta, MsDataScan dataScan, string dataFilePath,
            int precursorChargeState, List<Product> products)
        {
            List<Tuple<PeptideWithSetModifications, List<MatchedFragmentIon>, List<Product>>> variantsTuplesList = new();
            
            foreach (var peptide in peptideVariants)
            {
                peptide.Fragment(dataScan.DissociationType ?? DissociationType.HCD, FragmentationTerminus.Both, products);

                var match = MetaMorpheusEngine.MatchFragmentIons(new Ms2ScanWithSpecificMass(dataScan,
                    dataScan.SelectedIonMZ.Value,
                    precursorChargeState, dataFilePath, new CommonParameters()), products, new CommonParameters());

                variantsTuplesList.Add(
                    new Tuple<PeptideWithSetModifications, List<MatchedFragmentIon>, List<Product>>(
                        peptide, match, products));
            }

            var groupedSortedPeptides = variantsTuplesList.GroupBy(x => x.Item1.AllModsOneIsNterminus);

            List<Tuple<int[], int[]>> BYMatchesFromEachGroup = new();

            foreach (var group in groupedSortedPeptides)
            {
                var (encodedB, encodedY) =
                    FragmentMatchEncoder(group.First().Item1, group.First().Item2, group.First().Item3);
                BYMatchesFromEachGroup.Add(new Tuple<int[], int[]>(encodedB, encodedY));
            }

            return new Tuple<PeptideWithSetModifications, List<MatchedFragmentIon>, double, int>(null, null, 1, 2);
        }

        /// <summary>
        /// Decide either to remove combos with that amino acid position or stay with it to reduce unnecessary computation time. 
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="matchedFragmentIons"></param>
        /// <param name="inSilicoProducts"></param>
        /// <returns></returns>
        public (int[], int[]) FragmentMatchEncoder(PeptideWithSetModifications peptide, List<MatchedFragmentIon> matchedFragmentIons, List<Product> inSilicoProducts)
        {
            //Get arrays of ion fragments to compare
            var temp = matchedFragmentIons.Where(ion => ion.NeutralTheoreticalProduct.ProductType == ProductType.y)
                .Select(x => x.NeutralTheoreticalProduct.ProductType.ToString() + x.NeutralTheoreticalProduct.FragmentNumber);
            var matchedYFragmentIons = matchedFragmentIons.Where(ion => ion.NeutralTheoreticalProduct.ProductType == ProductType.y)
                .Select(x => x.NeutralTheoreticalProduct.ProductType.ToString() + x.NeutralTheoreticalProduct.FragmentNumber);
            var matchedBFragmentIons = matchedFragmentIons.Where(ion => ion.NeutralTheoreticalProduct.ProductType == ProductType.b)
                .Select(x => x.NeutralTheoreticalProduct.ProductType.ToString() + x.NeutralTheoreticalProduct.FragmentNumber);

            var inSilicoYFragmentIons = inSilicoProducts.Where(ion => ion.ProductType == ProductType.y)
                .Select(x => x.ProductType.ToString() + x.FragmentNumber);
            var inSilicoBFragmentIons = inSilicoProducts.Where(ion => ion.ProductType == ProductType.b)
                .Select(x => x.ProductType.ToString() + x.FragmentNumber);

            var modLocations = peptide.AllModsOneIsNterminus; //get the mod and location in peptide

            int[] bFragmentEncodedMatch = FormHotOneEncodedMatchArray(matchedBFragmentIons.ToArray(), inSilicoBFragmentIons.ToArray());
            int[] yFragmentEncodedMatch = FormHotOneEncodedMatchArray(matchedYFragmentIons.ToArray(), inSilicoYFragmentIons.ToArray());

            return (bFragmentEncodedMatch, yFragmentEncodedMatch);
        }

        /// <summary>
        /// Get hits of b and y ions
        /// </summary>
        /// <param name="bFragmentEncodedMatch"></param>
        /// <param name="yFragmentEncodedMatch"></param>
        /// <returns></returns>
        private (int, int) GetBYScores(int[] bFragmentEncodedMatch, int[] yFragmentEncodedMatch)
        {
            int bScore = 0;
            int yScore = 0;

            foreach (var b in bFragmentEncodedMatch)
            {
                if (b == 1) bScore++;
            }

            foreach (var y in yFragmentEncodedMatch)
            {
                if (y == 1) yScore++;
            }

            return (bScore, yScore);
        }

        /// <summary>
        /// array representing all positions in products, 1 is matched, 0 is not matched
        /// </summary>
        /// <returns></returns>
        private int[] FormHotOneEncodedMatchArray(string[] matchedArray, string[] productArray)
        {
            int[] zerosArray = new int[productArray.Count()];

            for (int i = 0; i < productArray.Count(); i++)
            {
                if (i > matchedArray.Count() - 1)
                    break;
                
                if (productArray[i] == matchedArray[i])
                    zerosArray[i] = 1;
            }
            
            return zerosArray;
        }

        /// <summary>
        /// Compares both matched and insilico ion products and returns the index of the last matching fragment present in both. [b,y respectively]
        /// </summary>
        /// <param name="inSilicoProductsArray"></param>
        /// <param name="matchedProductsArray"></param>
        private int CompareIonMatches(string[] inSilicoProductsArray, string[] matchedProductsArray)
        {
            
            for (int i = 0; i < matchedProductsArray.Length; i++)
            {
                if (matchedProductsArray[i] != inSilicoProductsArray[i])
                {
                    return i + 1;
                }
            }

            //returns -1 if all fragment ions matched!
            return -1;
        }

        private void SaveAsJson(string savingJsonPath)
        {
            List<MultipleSearchResults> listOfResults = new();

            foreach (var peptide in peptideGroup)
            {
                listOfResults.Add(new MultipleSearchResults()
                {
                    AccessionNumber = peptide.Item1.Protein.Accession,
                    BaseSequence = peptide.Item1.BaseSequence,
                    FullSequence = peptide.Item1.FullSequence,
                    IonMatchedCount = peptide.Item2.Count,
                    IsDecoy = peptide.Item1.Protein.IsDecoy,
                    MassErrorDa = peptide.Item2.Select(x => x.MassErrorDa).ToArray(),
                    MassErrorPpm = peptide.Item2.Select(x => x.MassErrorPpm).ToArray(),
                    MatchedIonCharge = peptide.Item2.Select(x => x.Charge).ToArray(),
                    MatchedIons = peptide.Item2.Select(x => x.Annotation).ToArray(),
                    MatchedMz = peptide.Item2.Select(x => x.Mz).ToArray(),
                    MonoisotopicMass = peptide.Item1.MonoisotopicMass,
                    PeptideLength = peptide.Item1.Length,
                    MostAbundantMonoisotopicMass = peptide.Item1.MostAbundantMonoisotopicMass,
                    TheoricalMz = peptide.Item2.Select(x => x.NeutralTheoreticalProduct.NeutralMass).ToArray(),
                    Modifications = String.Join(", ", peptide.Item1.AllModsOneIsNterminus
                        .Select(x => $"{x.Value.IdWithMotif}").ToArray()),
                    SequenceCoverage = peptide.Item2.Count,
                    MMmatch = peptide.Item3,
                    Charge = peptide.Item4
                });
            }
            var options = new JsonSerializerOptions { WriteIndented = true };
            string jsonString = JsonSerializer.Serialize(listOfResults.Select(x => x), options);
            File.WriteAllText($@"{savingJsonPath}", jsonString);
        }
        public static (List<DataTable>, List<MultipleSearchResults>) Run(MultipleSearchEngine engine,
            List<FilteredPsmTSV> psmList, List<Modification> fixedMods,
            MsDataFile dataFile, int maxNumOfMods, string pathToSaveResults)
        {


            List<KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>> tempMatches = new();


            Parallel.ForEach(psmList, psm =>
            {
                IEnumerable<PeptideWithSetModifications> peptideAsProtein = MultipleSearchEngine.GetPeptideAsProtein(psm, fixedMods);
                var matches = engine.GetPeptideFragmentIonsMatches(psm, dataFile, fixedMods);
                if (matches == null)
                    return;
                tempMatches.Add(matches);
            });

            var results = tempMatches.OrderByDescending(x => x.Value.Count)
                .GroupBy(x => x.Key.BaseSequence).ToList();


            List<DataTable> proteinGroupsTables = new();
            List<MultipleSearchResults> searchResults = new();
            Parallel.ForEach(results, result =>
            {
                var table = new DataTable();
                table.TableName = result.Key;
                foreach (var feature in typeof(MultipleSearchResults).GetProperties())
                {
                    table.Columns.Add(new DataColumn(feature.Name));
                }
                foreach (var peptide in result)
                {
                    var individualPeptide = new MultipleSearchResults()
                    {
                        AccessionNumber = peptide.Key.Protein.Accession,
                        BaseSequence = peptide.Key.Protein.BaseSequence,
                        FullSequence = peptide.Key.FullSequence,
                        IonMatchedCount = peptide.Value.Count,
                        Modifications = String.Join(", ", peptide.Key.AllModsOneIsNterminus
                            .Select(x => $"{x.Value.IdWithMotif}").ToArray()),
                        SequenceCoverage = peptide.Value.Count > 0 ? (peptide.Key.BaseSequence.Length * 2) / peptide.Value.Count : 0,
                        IsDecoy = peptide.Key.Protein.IsDecoy,
                        MassErrorDa = peptide.Value.Select(x => x.MassErrorDa).ToArray(),
                        MassErrorPpm = peptide.Value.Select(x => x.MassErrorPpm).ToArray(),
                        MatchedIonCharge = peptide.Value.Select(x => x.Charge).ToArray(),
                        MatchedIons = peptide.Value.Select(x => x.Annotation).ToArray(),
                        MatchedMz = peptide.Value.Select(x => x.Mz).ToArray(),
                        TheoricalMz = peptide.Value.Select(x => x.NeutralTheoreticalProduct.NeutralMass).ToArray(),
                        MonoisotopicMass = peptide.Key.MonoisotopicMass,
                        MostAbundantMonoisotopicMass = peptide.Key.MostAbundantMonoisotopicMass,
                        PeptideLength = peptide.Key.Protein.Length
                    };
                    DataRow row = table.NewRow();
                    row[0] = individualPeptide.BaseSequence;
                    row[1] = individualPeptide.SequenceCoverage;
                    row[2] = individualPeptide.IonMatchedCount;
                    row[3] = individualPeptide.Modifications;
                    row[4] = individualPeptide.FullSequence;
                    row[5] = individualPeptide.AccessionNumber;
                    row[6] = individualPeptide.PeptideLength;
                    row[7] = individualPeptide.MonoisotopicMass;
                    row[8] = individualPeptide.MostAbundantMonoisotopicMass;
                    row[9] = individualPeptide.IsDecoy;
                    row[10] = String.Join(", ", individualPeptide.MatchedIons);
                    row[11] = String.Join(", ", individualPeptide.MatchedIonCharge);
                    row[12] = String.Join(", ", individualPeptide.TheoricalMz);
                    row[13] = String.Join(", ", individualPeptide.MatchedMz);
                    row[14] = String.Join(", ", individualPeptide.MassErrorPpm);
                    row[15] = String.Join(", ", individualPeptide.MassErrorDa);

                    table.Rows.Add(row);
                    searchResults.Add(individualPeptide);

                }
                proteinGroupsTables.Add(table);
            });

            var options = new JsonSerializerOptions { WriteIndented = true };
            string jsonString = JsonSerializer.Serialize(searchResults.Select(x => x), options);
            File.WriteAllText($@"{pathToSaveResults}", jsonString);

            return (proteinGroupsTables, searchResults);
            //proteinGroups.ItemsSource = proteinGroupsTables.AsEnumerable();
        }

        private static KeyValuePair<double, Modification[]>[] SortAndEliminateDuplicates(
            KeyValuePair<double, Modification[]>[] combinationsFromDatabase)
        {
            Parallel.ForEach(combinationsFromDatabase, combo =>
            {
                combo = new KeyValuePair<double, Modification[]>(combo.Key,
                    combo.Value.OrderBy(mod => mod.IdWithMotif).ToArray());
            });
            return combinationsFromDatabase.Distinct().OrderBy(x => x.Key).ToArray();
        }

        public List<KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>> GetPeptideFragmentIonsMatches(
            FilteredPsmTSV psm, MsDataFile dataFile, List<Modification> fixedMods)
        {
            List<KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>> results = new();

            var spectrum = dataFile.GetOneBasedScan(int.Parse(psm.ScanNumber));
            //var neutralMass = spectrum.SelectedIonMZ.Value.ToMass(spectrum.SelectedIonChargeStateGuess.Value);

            var peptide = GetPeptideWithMods(psm, fixedMods);

            var deltaMass = GetDeltaMass(psm, peptide.First());

            //if (deltaMass < 1)
            //    return results;

            var probableMods =
                GetCombinationsThatFitDelta(deltaMass);

            //var products = new List<Product>();
            int id = 0;
            foreach (var mod in probableMods)
            {
                var products = new List<Product>();
                var ptm = GetPeptideWithMods(psm, fixedMods, mod.ToList());
                foreach (var variant in ptm)
                {
                    variant.Fragment(spectrum.DissociationType ?? DissociationType.HCD, FragmentationTerminus.Both,
                        products);

                    var match = MetaMorpheusEngine.MatchFragmentIons(new Ms2ScanWithSpecificMass(spectrum,
                        double.Parse(psm.PrecursorMass),
                        spectrum.SelectedIonChargeStateGuess.Value, dataFile.FilePath,
                        new CommonParameters()), products, new CommonParameters());
                    if (id == 0)
                    {

                        results.Add(new KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>(variant, match));
                        id = id + 1;
                    }
                    else
                    {
                        if (match.Count < results[results.Count - 1].Value.Count)
                        {
                            results.Add(new KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>(variant, match));
                            id = id + 1;
                        }
                    }

                    if (GetCombinationsThatFitDelta(GetDeltaMass(psm, variant)).Count > 0)
                    {

                    }
                }
            }
            return results;
        }


        /// <summary>
        /// Returns the mass difference between chosen precursor ion and given peptide.
        /// </summary>
        /// <returns></returns>
        private static double GetDeltaMass(FilteredPsmTSV psm, PeptideWithSetModifications peptide)
        {
            //var spectrum = dataFile.GetOneBasedScan(int.Parse(psm.ScanNumber));
            var precursorMass = double.Parse(psm.PrecursorMass);
            //var precursorIonMass = spectrum.SelectedIonMonoisotopicGuessMz.Value * spectrum.SelectedIonChargeStateGuess.Value;

            return precursorMass - peptide.MonoisotopicMass;

        }

        /// <summary>
        /// Returns the Combination of mods that are possible. Uses binary search on the IOrderedEnumerable of mods.
        /// </summary>
        /// <param name="possibleMods"></param>
        /// <param name="peptide"></param>
        /// <param name="psm"></param>
        /// <param name="dataFile"></param>
        /// <returns></returns>
        public List<List<Modification>> GetCombinationsThatFitDelta(double deltaMass)
        {
            var tolerance = new PpmTolerance(50);

            //var massArray = CombinationsFromDatabase.OrderBy(x => x.Key).Select(x => x.Key).ToArray();

            var temp1 = tolerance.GetMaximumValue(deltaMass);
            var temp2 = tolerance.GetMinimumValue(deltaMass);

            var maxIndex = ~Array.BinarySearch(MassArray, tolerance.GetMaximumValue(deltaMass));
            var minIndex = ~Array.BinarySearch(MassArray, tolerance.GetMinimumValue(deltaMass));

            if (maxIndex == minIndex && maxIndex >= CombinationOfModifications.Count || minIndex >= CombinationOfModifications.Count)
                return new List<List<Modification>>();

            int[] rangeIndex = Enumerable.Range(minIndex, Math.Abs(maxIndex - minIndex)).ToArray();

            if (rangeIndex.Length == 0)
                rangeIndex = new[] { maxIndex };

            var rangeOfPossibleMods = rangeIndex.Select(x => CombinationOfModifications[x]);

            return rangeOfPossibleMods.ToList();
        }

        /// <summary>
        /// Returns collection of custom psms. These are part of a custom class for development, not for production.
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public static List<FilteredPsmTSV> ReadFilteredPsmTSVShort(string path)
        {
            List<FilteredPsmTSV> filteredList = new List<FilteredPsmTSV>();

            using (var reader = new StreamReader(path))
            {
                reader.ReadLine();
                string[] lineCheck;
                while (reader.EndOfStream == false)
                {
                    var line = reader.ReadLine().Split('\t');
                    FilteredPsmTSV filteredPsm = new FilteredPsmTSV(line);
                    filteredList.Add(filteredPsm);
                }
            }

            return filteredList;
        }

        private static IEnumerable<PeptideWithSetModifications> GetPeptideWithMods(FilteredPsmTSV psm, List<Modification> fixedMods = null, List<Modification> variableMods = null)
        {
            var protein = new Protein(psm.BaseSeq, psm.ProteinAccession);
            return protein.Digest(new DigestionParams("top-down"), fixedMods, variableMods);
        }

        /// <summary>
        /// Returns peptide threated as a protein using top-down digestion. Used for GPTMD sequence lead.
        /// </summary>
        /// <param name="psm"></param>
        /// <param name="fixedMods"></param>
        /// <returns></returns>
        public static IEnumerable<PeptideWithSetModifications> GetPeptideAsProtein(FilteredPsmTSV psm, List<Modification> fixedMods = null)
        {
            var peptideProteinDigest =
                new Protein(psm.BaseSeq, psm.ProteinAccession).Digest(new DigestionParams(protease: "top-down"),
                    allKnownFixedModifications: fixedMods,
                    variableModifications: new List<Modification>());

            return peptideProteinDigest;
        }

        /// <summary>
        /// Returns mods from a toml containing common modifications. 
        /// </summary>
        /// <param name="gptmdToml"></param>
        /// <returns></returns>
        public static IEnumerable<Modification> GetModsFromGptmdThing(string gptmdToml = @"Task1-GPTMDTaskconfig.toml")
        {
            var task = Toml.ReadFile<GptmdTask>(gptmdToml,
                MetaMorpheusTask.tomlConfig);

            var mods = GlobalVariables.AllModsKnownDictionary;

            foreach (var (item1, item2) in task.GptmdParameters.ListOfModsGptmd)
            {
                if (mods.TryGetValue(item2, out Modification mod))
                {
                    yield return mod;
                }
            }
        }
        public static void Mc(List<Modification> sortedListOfModsToAdd, ref List<List<Modification>> listOfModCombinations, int maxNumberOfModsInGroup, bool   allModsFromOneToN)
        {
            List<List<Modification>> newModLists = new();
            if (maxNumberOfModsInGroup == 0)
            {
                foreach (var modificationToAdd in sortedListOfModsToAdd)
                {
                    newModLists.Add(new List<Modification>() { modificationToAdd });
                }
                newModLists = newModLists.DistinctBy(n => ModListNameString(n)).ToList();
                listOfModCombinations = newModLists;
            }
            else
            {
                Mc(sortedListOfModsToAdd, ref listOfModCombinations, maxNumberOfModsInGroup - 1, allModsFromOneToN);
                newModLists.Clear();
                foreach (var modList in listOfModCombinations.Where(c => c.Count == (maxNumberOfModsInGroup - 1)))
                {
                    foreach (var modificationToAdd in sortedListOfModsToAdd)
                    {
                        List<Modification> newModList = modList.ToList();
                        newModList.Add(modificationToAdd);
                        newModLists.Add(newModList.OrderBy(n => n.IdWithMotif).ToList());
                    }
                }
                newModLists = newModLists.DistinctBy(n => ModListNameString(n)).ToList();
                listOfModCombinations.AddRange(newModLists);
                if (!allModsFromOneToN)
                {
                    listOfModCombinations = listOfModCombinations.Where(c => c.Count == maxNumberOfModsInGroup).ToList();
                }
            }

        }
        public static string ModListNameString(List<Modification> list)
        {
            return String.Join("", list.Select(n => n.IdWithMotif));
        }
    }
}
