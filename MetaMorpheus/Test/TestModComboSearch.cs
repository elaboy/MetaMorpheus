using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using ExCSS;
using iText.Kernel.Pdf.Canvas.Parser.ClipperLib;
using iText.StyledXmlParser.Jsoup.Helper;
using TaskLayer;
using ThermoFisher.CommonCore.Data;
using UsefulProteomicsDatabases;

namespace Test
{
    public sealed class TestModComboSearch
    {
        private const string _filteredPsm = @"D:\08-30-22_bottomup\example.psmtsv";

        //[Test]
        //public void TestMultiModSearch()
        //{
        //    var mods =
        //        Loaders.LoadUnimod(
        //                @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
        //            .ToList();

        //    int maxToAdd = 7;
        //    var fixedMods = new List<Modification>();
        //    fixedMods.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")));

        //    var commonVariableMods = MultiModSearch.GetModsFromGptmdThing().ToList();

        //    var psmList =
        //        MultiModSearch.ReadFilteredPsmTSVShort(_filteredPsm);

        //    var msDataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML").GetAllScansList()
        //        .Where(x => x.MsnOrder == 2).ToArray();
        //    var dataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML");
        //    List<KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>> allMatches = new();
        //    List<FilteredPsmTSV> notBetter = new();

        //    for (int i = 1; i < maxToAdd + 1; i++)
        //    {
        //        var modSearch = new MultiModSearch(commonVariableMods, i);

        //        List<KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>> tempMatches = new();

        //        Parallel.ForEach(psmList, psm =>
        //        {
        //            var peptideAsProtein = MultiModSearch.GetPeptideAsProtein(psm, fixedMods);
        //            var matches = modSearch.GetPeptideFragmentIonsMatches(psm, dataFile, fixedMods);
        //            var top = matches.OrderByDescending(x => x.Value.Count);
        //            if (matches.Count is 0)
        //                return;
        //            var selected = top.First();

        //            tempMatches.Add(selected);
        //        });

        //        var orderedMatches = tempMatches.OrderByDescending(x => x.Value.Count())
        //            .GroupBy(x => x.Key.Protein.BaseSequence)
        //            .Select(x => new Comparing()
        //            {
        //                mmMatchCount = int.Parse(psmList
        //                    .Find(psm => x.Key
        //                        .Equals(psm.BaseSeq)).MatchedIonCounts),
        //                myMatchCount = x.OrderByDescending(x => x.Value.Count)
        //                    .First().Value.Count,
        //                peptideSequence = x.Key
        //            });

        //        foreach (var match in orderedMatches)
        //        {
        //            if(match.myMatchCount < int.Parse(psmList.Find(psm => psm.BaseSeq.Equals(match.peptideSequence)).MatchedIonCounts))
        //            {
        //                notBetter.Add(psmList.Find(psm => psm.BaseSeq.Equals(match.peptideSequence)));
        //            }
        //        }

        //        psmList = notBetter;

        //        tempMatches.ForEach(x => allMatches
        //                .Add(new KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>(x.Key, x.Value)));
        //    }

        //    int zero = 0;
        //}

        [Test]
        public void TestModCombinations()
        {
            var commonVariableMods = MultiModSearch.GetModsFromGptmdThing().ToList();


            List<Modification> chosenSubsetOfMods = commonVariableMods.ToArray()[0..5].ToList();
            List<List<Modification>> completeListOfModCombinations = new();
            Mc(chosenSubsetOfMods, ref completeListOfModCombinations, 3, false);
            File.WriteAllLines(@"C:\Users\Edwin\Desktop\mods.txt", completeListOfModCombinations.Select(n => ModListNameString(n)).ToArray());
            Assert.AreEqual(1, 0);

            Assert.IsTrue(false);
        }

        [Test]
        public void TestModCombinationsMatch()
        {

            var mods =
                Loaders.LoadUnimod(
                        @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                    .ToList();
            int maxToAdd = 5;

            var fixedMods = new List<Modification>();
            fixedMods.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")));


            var commonVariableMods = MultiModSearch.GetModsFromGptmdThing().ToList();


            List<Modification> commonBiologycalMods = GlobalVariables.AllModsKnown.OfType<Modification>()
                .Where(mod => mod.ModificationType.Contains("Common Biological")).ToList();

            var engine = new MultiModSearch(commonBiologycalMods, maxToAdd, true);

            var psmList =
                MultiModSearch.ReadFilteredPsmTSVShort(_filteredPsm);

            var msDataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML").GetAllScansList()
                .Where(x => x.MsnOrder == 2).ToArray();
            var dataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML");


            List<KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>> tempMatches = new();

            Parallel.ForEach(psmList, psm =>
            {
                var peptideAsProtein = MultiModSearch.GetPeptideAsProtein(psm, fixedMods);
                var matches = engine.GetPeptideFragmentIonsMatches(psm, dataFile, fixedMods);

                tempMatches.Add(matches);
            });

            var results = tempMatches.OrderByDescending(x => x.Value.Count)
                .GroupBy(x => x.Key.BaseSequence);

            List<List<MultiModSearchResults>> peptideGroups = new();

            Parallel.ForEach(results, result =>
            {
                List<MultiModSearchResults> groupedPeptides = new();
                foreach (var peptide in result)
                {
                    groupedPeptides.Add(new MultiModSearchResults()
                    {
                        AccessionNumber = peptide.Key.Protein.Accession,
                        BaseSequece = peptide.Key.Protein.BaseSequence,
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
                    });
                }
                peptideGroups.Add(groupedPeptides);
            });
            var options = new JsonSerializerOptions { WriteIndented = true };
            string jsonString = JsonSerializer.Serialize(peptideGroups.Select(x => x), options);
            File.WriteAllText(@"C:\Users\Edwin\Desktop\jsonExample", jsonString);

        }
        /// <summary>
        /// Groups of base sequences will be stored in this objecto to serialize into a xml database
        /// </summary>
        public class MultiModSearchResults
        {
            public string BaseSequece { get; set; }
            public string AccessionNumber { get; set; }
            public int PeptideLength { get; set; }
            public double MonoisotopicMass { get; set; }
            public double MostAbundantMonoisotopicMass { get; set; }
            public bool IsDecoy { get; set; }
            public string[] MatchedIons { get; set; }
            public int[] MatchedIonCharge { get; set; }
            public double[] TheoricalMz { get; set; }
            public double[] MatchedMz { get; set; }
            public double[] MassErrorPpm { get; set; }
            public double[] MassErrorDa { get; set; }
        }

        [Test]
        public static void Bubba()
        {
            List<Modification> allAvailableMods = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "CommonBiological.txt"), out var errors).ToList();
            
            List<Modification> chosenSubsetOfMods = allAvailableMods.ToArray()[0..5].ToList();
            List<List<Modification>> completeListOfModCombinations = new();
            Mc(chosenSubsetOfMods, ref completeListOfModCombinations, 3, true);
            File.WriteAllLines(@"E:\junk\mods.txt", completeListOfModCombinations.Select(n => ModListNameString(n)).ToArray());
            Assert.AreEqual(1, 0);
        }
        public static void Mc(List<Modification> sortedListOfModsToAdd, ref List<List<Modification>> listOfModCombinations, int maxNumberOfModsInGroup, bool allModsFromOneToN)
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


    public class Comparing
    {
        public string peptideSequence { get; set; }
        public int mmMatchCount {get; set; }
        public int myMatchCount { get; set; }

    }

    public class MultiModSearch
    {
        public  KeyValuePair<double, Modification[]>[] CombinationsFromDatabase { get; set; }
        public  IOrderedEnumerable<IGrouping<string, Modification>> ModificationGroups { get; set; }
        public double[] MassArray { get; set; }

        public List<List<Modification>> CombinationOfModifications { get; set; }
        public List<KeyValuePair<double, Modification[]>> CombinationsWithAddedMass { get; set; }

        //public MultiModSearch(List<Modification> listOfMods, int numberOfVariableMods)
        //{
        //    ModificationGroups = SetModsGroup(listOfMods);
        //    CombinationsFromDatabase = SetKCombinations(ModificationGroups, numberOfVariableMods).ToArray();
        //    CombinationsFromDatabase = SortAndEliminateDuplicates(CombinationsFromDatabase);
        //    MassArray = CombinationsFromDatabase.Select(x => x.Key).ToArray();
        //}

        public MultiModSearch(List<Modification> listOfMods, int numberOfVariableMods, bool allCombos)
        {
            List<List<Modification>> comboList = new();
            CombinationsWithAddedMass = new();
            Mc(listOfMods, ref comboList, numberOfVariableMods, allCombos);
            CombinationOfModifications = comboList;
            CombinationsWithAddedMass.Add(CombinationOfModifications.Select(x => new KeyValuePair<double, Modification[]>(
                key: x.Select(x => x.MonoisotopicMass.Value).Sum(), value: x.Select(x => x).ToArray())));
            CombinationsWithAddedMass = CombinationsWithAddedMass.OrderBy(x => x.Key).ToList();
            MassArray = CombinationsWithAddedMass.Select(x => x.Key).ToArray();
        }


        private static KeyValuePair<double, Modification[]>[] SortAndEliminateDuplicates(KeyValuePair<double, Modification[]>[] combinationsFromDatabase)
        {
            Parallel.ForEach(combinationsFromDatabase, combo =>
            {
                combo = new KeyValuePair<double, Modification[]>(combo.Key,
                    combo.Value.OrderBy(mod => mod.IdWithMotif).ToArray());
            });
            return combinationsFromDatabase.Distinct().OrderBy(x => x.Key).ToArray();
        }

        public List<KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>> GetPeptideFragmentIonsMatches(FilteredPsmTSV psm, MsDataFile dataFile, List<Modification> fixedMods)
        {
            List<KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>> results = new();

            var spectrum = dataFile.GetOneBasedScan(int.Parse(psm.ScanNumber));
            //var neutralMass = spectrum.SelectedIonMZ.Value.ToMass(spectrum.SelectedIonChargeStateGuess.Value);

            var peptide = GetPeptideWithMods(psm, fixedMods);

            var deltaMass = GetDeltaMass(psm, peptide.First());

            if (deltaMass < 1)
                return results;

            var probableMods =
                GetCombinationsThatFitDelta(deltaMass);

            //var products = new List<Product>();
            int id = 0;
            foreach(var mod in probableMods)
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

                    results.Add(new KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>(variant, match));
                    id = id + 1;
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
        private List<List<Modification>> GetCombinationsThatFitDelta(double deltaMass)
        {
            var tolerance = new PpmTolerance(15);

            //var massArray = CombinationsFromDatabase.OrderBy(x => x.Key).Select(x => x.Key).ToArray();

            var temp1 = tolerance.GetMaximumValue(deltaMass);
            var temp2 = tolerance.GetMinimumValue(deltaMass);

            var maxIndex = ~Array.BinarySearch(MassArray, tolerance.GetMaximumValue(deltaMass));
            var minIndex = ~Array.BinarySearch(MassArray, tolerance.GetMinimumValue(deltaMass));

            if(maxIndex == minIndex && maxIndex >= CombinationOfModifications.Count || minIndex >= CombinationOfModifications.Count)
                return new List<List<Modification>>();

            int[] rangeIndex = Enumerable.Range(minIndex, Math.Abs(maxIndex - minIndex)).ToArray();

            if (rangeIndex.Length == 0)
                rangeIndex = new[] { maxIndex};

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

        private static IEnumerable<KeyValuePair<double, Modification[]>> SetKCombinations(
            IEnumerable<IGrouping<string, Modification>> mods, int numberOfMods)
        {
            var combination = new List<KeyValuePair<double, Modification[]>>();

            switch (numberOfMods)
            {
                case 1:
                    var oneMod = from mod in mods.AsParallel()
                                 select new KeyValuePair<double, Modification[]>(
                                     mod.First().MonoisotopicMass.Value, new[] { mod.First() });
                    return oneMod;

                case 2:
                    var twoMods = from mod in mods.AsParallel()
                                  from mod2 in mods.AsParallel()
                                  select new KeyValuePair<double, Modification[]>(
                                      mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value,
                                      new[] { mod.First(), mod2.First() }.OrderBy(x => x.IdWithMotif).ToArray());
                    return twoMods;

                case 3:
                    var threeMods = from mod in mods.AsParallel()
                                    from mod2 in mods.AsParallel()
                                    from mod3 in mods.AsParallel()
                                    select new KeyValuePair<double, Modification[]>(
                                        mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value + mod3.First().MonoisotopicMass.Value,
                                        new[] { mod.First(), mod2.First(), mod3.First() }.OrderBy(x => x.IdWithMotif).ToArray());
                    return threeMods; 

                case 4:
                    var fourMods = from mod in mods.AsParallel()
                                   from mod2 in mods.AsParallel()
                                   from mod3 in mods.AsParallel()
                                   from mod4 in mods.AsParallel()
                                   select new KeyValuePair<double, Modification[]>(
                                       mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value + mod3.First().MonoisotopicMass.Value
                                       + mod4.First().MonoisotopicMass.Value,
                                       new[] { mod.First(), mod2.First(), mod3.First(), mod4.First() }.OrderBy(x => x.IdWithMotif).ToArray());
                    return fourMods;

                case 5:
                    var fiveMods = from mod in mods.AsParallel()
                                   from mod2 in mods.AsParallel()
                                   from mod3 in mods.AsParallel()
                                   from mod4 in mods.AsParallel()
                                   from mod5 in mods.AsParallel()
                                   select new KeyValuePair<double, Modification[]>(
                                       mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value + mod3.First().MonoisotopicMass.Value
                                       + mod4.First().MonoisotopicMass.Value + mod5.First().MonoisotopicMass.Value,
                                       new[] { mod.First(), mod2.First(), mod3.First(), mod4.First(), mod5.First() }.OrderBy(x => x.IdWithMotif).ToArray());
                    return fiveMods;

                case 6:
                    var sixMods = from mod in mods.AsParallel()
                                  from mod2 in mods.AsParallel()
                                  from mod3 in mods.AsParallel()
                                  from mod4 in mods.AsParallel()
                                  from mod5 in mods.AsParallel()
                                  from mod6 in mods.AsParallel()
                                  select new KeyValuePair<double, Modification[]>(
                                      mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value + mod3.First().MonoisotopicMass.Value
                                      + mod4.First().MonoisotopicMass.Value + mod5.First().MonoisotopicMass.Value + mod6.First().MonoisotopicMass.Value,
                                      new[] { mod.First(), mod2.First(), mod3.First(), mod4.First(), mod5.First(), mod6.First() }.OrderBy(x => x.IdWithMotif).ToArray());
                    return sixMods;

                case 7:
                    var sevenMods = from mod in mods.AsParallel()
                                    from mod2 in mods.AsParallel()
                                    from mod3 in mods.AsParallel()
                                    from mod4 in mods.AsParallel()
                                    from mod5 in mods.AsParallel()
                                    from mod6 in mods.AsParallel()
                                    from mod7 in mods.AsParallel()
                                    select new KeyValuePair<double, Modification[]>(
                                        mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value + mod3.First().MonoisotopicMass.Value
                                        + mod4.First().MonoisotopicMass.Value + mod5.First().MonoisotopicMass.Value + mod6.First().MonoisotopicMass.Value
                                        + mod7.First().MonoisotopicMass.Value,
                                        new[] { mod.First(), mod2.First(), mod3.First(), mod4.First(), mod5.First(), mod6.First(), mod7.First() }.OrderBy(x => x.IdWithMotif).ToArray());
                    return sevenMods;
            }
            
            return combination.OrderBy(x => x.Key);
        }

        /// <summary>
        /// Sets the database and groups it by Modification ID. For example: Acetylations.
        /// </summary>
        /// <param name="listOfMods"></param>
        /// <returns></returns>
        private static IOrderedEnumerable<IGrouping<string, Modification>> SetModsGroup(List<Modification> listOfMods)
        {
            var commonModsFromToml = GetModsFromGptmdThing().ToList();

            var mods =
                Loaders.LoadUnimod(
                        @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                    .ToList();

            commonModsFromToml.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C"))); //adds this mod into the toml database (it's a fixed mod)

            return commonModsFromToml.GroupBy(modId => modId.IdWithMotif).OrderBy(x => x.Key);
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
        public static void Mc(List<Modification> sortedListOfModsToAdd, ref List<List<Modification>> listOfModCombinations, int maxNumberOfModsInGroup, bool allModsFromOneToN)
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
