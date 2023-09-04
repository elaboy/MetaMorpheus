using MassSpectrometry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using Easy.Common.Extensions;
using System.Threading.Tasks;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using MzLibUtil;
using System.IO;
using FlashLFQ;
using Microsoft.ML.Transforms;

namespace EngineLayer.CombinatorialSearch
{
    public class CSEngine : MetaMorpheusEngine
    {

        public List<List<Modification>> CombinationOfModifications { get; set; }
        public List<KeyValuePair<double, Modification[]>> CombinationsWithAddedMass { get; set; }
        public double[] MassArray { get; set; }
        public List<Tuple<double, Protein, MsDataScan, double, int>> ProteinListInferedFromGPTMD { get; private set; }
        private List<Dictionary<int, Modification>> ModsToIgnore { get; set; }
        private List<Modification> FixedMods { get; set; }
        private List<FilteredPsmTSV> PsmsList { get; set; }
        public List<PeptideSpectralMatch> Psms { get; set; }
        private MsDataFile MsDataFile { get; set; }

        /// <summary>
        /// Constructor for the engine.
        /// </summary>
        /// <param name="psmList"></param>
        /// <param name="listOfMods"></param>
        /// <param name="numberOfVariableMods"></param>
        /// <param name="fixedMods"></param>
        /// <param name="commonParameters"></param>
        /// <param name="fileSpecificParameters"></param>
        /// <param name="nestedIds"></param>
        public CSEngine(List<PeptideSpectralMatch> peptideSpectralMatches, List<List<Modification>> listOfModCombos,
            List<KeyValuePair<double, Modification[]>> numberOfVariableMods, List<Modification> fixedMods,
            CommonParameters commonParameters, MsDataFile dataFile, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            List<string> nestedIds) 
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            Psms = peptideSpectralMatches;
            CombinationOfModifications = listOfModCombos;
            CombinationsWithAddedMass = numberOfVariableMods;
            //List<List<Modification>> comboList = new();
            //CombinationsWithAddedMass = new();
            ////CombinationBuilder(listOfMods, ref comboList, numberOfVariableMods, allCombos);
            ////CombinationOfModifications = comboList;
            //CombinationsWithAddedMass.Add(CombinationOfModifications.Select(x => new KeyValuePair<double, Modification[]>(
            //    key: x.Select(x => x.MonoisotopicMass.Value).Sum(), value: x.Select(x => x).ToArray())));
            //CombinationsWithAddedMass = CombinationsWithAddedMass.OrderBy(x => x.Key).ToList();
            //MassArray = CombinationsWithAddedMass.Select(x => x.Key).ToArray();
            //SetProteinsInfered(psmList, dataFile);
            ModsToIgnore = new List<Dictionary<int, Modification>>();
            //PsmsList = psmList;
            MsDataFile = dataFile;
            FixedMods = fixedMods;
        }

        /// <summary>
        /// Loads the Proteins into the engine. 
        /// </summary>
        /// <param name="psmList"></param>
        /// <param name="dataFile"></param>
        private void SetProteinsInfered(List<FilteredPsmTSV> psmList, MsDataFile dataFile)
        {
            //Todo: Change this weird tuple structure to something better both in performance and easier to use
            ProteinListInferedFromGPTMD = new();
            foreach(var psm in psmList)
            {
                ProteinListInferedFromGPTMD.Add(new Tuple<double, Protein, MsDataScan, double, int>(
                    double.Parse(psm.PrecursorMass),
                    new Protein(psm.BaseSeq, psm.ProteinAccession, psm.OrganismName),
                    dataFile.GetOneBasedScan(int.Parse(psm.ScanNumber)),
                    double.Parse(psm.Score), int.Parse(psm.Charge)));
            }
        }
        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Dictionary<string, HashSet<Tuple<int, Modification>>> modsUsedDictionary = new();
            lock (modsUsedDictionary)
            {
                //Parallel.ForEach(Psms, psm =>
                //{
                //    var peptide = new PeptideWithSetModifications(new Protein(psm.))
                //}


                Parallel.ForEach(ProteinListInferedFromGPTMD, protein =>
                {
                    Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> resultsFromSearch = new();

                    var peptide = protein.Item2
                        .Digest(new DigestionParams("top-down", 0), FixedMods,
                            new List<Modification>());

                    var deltaMass = protein.Item1 - peptide.First().MonoisotopicMass;

                    var possibleMods = GetCombinationsThatFitDelta(deltaMass);

                    var psm = PsmsList.Find(x => x.BaseSeq.Equals(protein.Item2.BaseSequence));

                    var msDataScan = protein.Item3;//MsDataFile.GetOneBasedScan(int.Parse(psm.ScanNumber));

                    foreach (var mod in possibleMods)
                    {
                        var bestCandidate = SearchForMods(FixedMods, possibleMods,
                            protein.Item2.Digest(new DigestionParams("top-down", 0),
                                FixedMods, mod),
                           msDataScan, psm, MsDataFile);

                        if (resultsFromSearch.Count == 0)
                        {
                            resultsFromSearch.Add(bestCandidate);
                            continue;
                        }

                        if (bestCandidate.Count == 0)
                        {
                            continue;
                        }

                        if (!resultsFromSearch.ContainsKey(bestCandidate.Keys.First()))
                        {
                            resultsFromSearch.Add(bestCandidate);
                        }

                        if (modsUsedDictionary.ContainsKey(protein.Item2.BaseSequence))
                        {
                            foreach (var peptideWithMods in bestCandidate)
                            {
                                foreach (var modInDict in peptideWithMods.Key.AllModsOneIsNterminus)
                                {
                                    modsUsedDictionary[protein.Item2.BaseSequence].Add(new Tuple<int, Modification>(
                                        modInDict.Key,
                                        modInDict.Value));

                                }
                            }
                        }
                        else
                        {
                            modsUsedDictionary.Add(protein.Item2.BaseSequence, new HashSet<Tuple<int, Modification>>());
                        }
                    }
                });
            }

            return new CSResults(this, modsUsedDictionary, 
                CombinationOfModifications.SelectMany(x => x).ToList(),
                ProteinListInferedFromGPTMD.Select(x => x.Item2).ToList());
        }

        /// <summary>
        /// Performs the search and returns a dictionaty of 
        /// </summary>
        /// <param name="fixedMods"></param>
        /// <param name="possibleModsThatFitDeltaMass"></param>
        /// <param name="peptidesToSearch"></param>
        /// <param name="scanForPsm1"></param>
        /// <param name="psm1"></param>
        /// <param name="dataFile"></param>
        /// <returns></returns>
        public Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> SearchForMods(List<Modification> fixedMods,
            List<List<Modification>> possibleModsThatFitDeltaMass, IEnumerable<PeptideWithSetModifications> peptidesToSearch,
            MsDataScan scanForPsm1, FilteredPsmTSV psm1, MsDataFile dataFile)
        {
            List<PeptideWithSetModifications> filteredPeptides = new();
            if (ModsToIgnore.Count > 0)
            {
                filteredPeptides = FilterPeptidesWithModsToIgnore(peptidesToSearch);
            }
            else
            {
                filteredPeptides = peptidesToSearch.ToList();
            }
            Dictionary<PeptideWithSetModifications, Tuple<List<MatchedFragmentIon>, List<MatchedFragmentIon>>> tempMatchesFragmentIons = new(); //b and y fragments tuple

            var bProducts = new List<Product>();
            var yProducts = new List<Product>();

            foreach (var peptide in filteredPeptides)
            {
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.N, bProducts);
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.C, yProducts);

                var bMatch = MetaMorpheusEngine.MatchFragmentIons(
                    new Ms2ScanWithSpecificMass(scanForPsm1, scanForPsm1.SelectedIonMZ.Value,
                        int.Parse(psm1.Charge), dataFile.FilePath,
                        new CommonParameters()), bProducts, new CommonParameters());

                var yMatch = MetaMorpheusEngine.MatchFragmentIons(
                    new Ms2ScanWithSpecificMass(scanForPsm1, scanForPsm1.SelectedIonMZ.Value,
                        int.Parse(psm1.Charge), dataFile.FilePath,
                        new CommonParameters()), yProducts, new CommonParameters());

                if (!tempMatchesFragmentIons.ContainsKey(peptide))
                    tempMatchesFragmentIons.Add(peptide, new Tuple<List<MatchedFragmentIon>, List<MatchedFragmentIon>>(bMatch, yMatch));
            }

            //Get rid of zero matching on either side
            var tempMatchesFragmentIonsWithoutZeros =
                tempMatchesFragmentIons.Where(x => x.Value.Item1.Count > 0 && x.Value.Item2.Count > 0);

            if (tempMatchesFragmentIonsWithoutZeros.Count() == 0)
                return new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();

            var goodMatches = AddNonCandidatesToModsToIgnore(tempMatchesFragmentIonsWithoutZeros);

            var resultsBByDescending = goodMatches.OrderByDescending(x => x.Value.Item1.Count);
            var resultsYByDescending = goodMatches.OrderByDescending(x => x.Value.Item2.Count);


            var bestModsFromResultsBandY = resultsYByDescending.First().Key.AllModsOneIsNterminus;
            foreach (var combo in resultsYByDescending.First().Key.AllModsOneIsNterminus)
            {
                if (!bestModsFromResultsBandY.ContainsKey(combo.Key))
                    bestModsFromResultsBandY.Add(combo.Key, combo.Value);
            }

            return GetBestPeptideMatch(bestModsFromResultsBandY, resultsBByDescending.First().Key.Protein, fixedMods,
                scanForPsm1, psm1, dataFile);
        }
        /// <summary>
        /// Filters the peptides so redundant matching is avoided.
        /// </summary>
        /// <param name="peptideWithSetModifications"></param>
        /// <returns></returns>
        private List<PeptideWithSetModifications> FilterPeptidesWithModsToIgnore(IEnumerable<PeptideWithSetModifications> peptideWithSetModifications)
        {
            List<PeptideWithSetModifications> filteredPeptides = new();

            foreach (var peptide in peptideWithSetModifications)
            {
                if (!ModsToIgnore.Contains(peptide.AllModsOneIsNterminus))
                    filteredPeptides.Add(peptide);
            }

            return filteredPeptides;
        }

        private Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> GetBestPeptideMatch(Dictionary<int, Modification> bestModsFromBYMatching, Protein protein,
            List<Modification> fixedMods, MsDataScan scanForPsm1, FilteredPsmTSV psm1, MsDataFile dataFile)
        {
            var secondRunResults = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();

            var secondRunPeptide = new PeptideWithSetModifications(protein,
                new DigestionParams("top-down", 0), 1, protein.BaseSequence.Length,
                CleavageSpecificity.Full, "", 0, bestModsFromBYMatching, 0);

            var secondRunDigested = secondRunPeptide.Protein.Digest(new DigestionParams("top-down", 0), fixedMods,
                bestModsFromBYMatching.Values.ToList());

            var secondProductsRun = new List<Product>();
            foreach (var product in secondRunDigested.Distinct())
            {
                product.Fragment(DissociationType.HCD, FragmentationTerminus.Both, secondProductsRun);

                var secondRunMatch = MetaMorpheusEngine.MatchFragmentIons(
                    new Ms2ScanWithSpecificMass(scanForPsm1, scanForPsm1.SelectedIonMZ.Value,
                        int.Parse(psm1.Charge), dataFile.FilePath,
                        new CommonParameters()), secondProductsRun, new CommonParameters());

                secondRunResults.Add(product, secondRunMatch);
            }

            var bestMatch = secondRunResults.OrderByDescending(x =>
                x.Value.Count).ThenByDescending(x => x.Value[0].NeutralTheoreticalProduct.FragmentNumber).First();

            var bestMatchDict = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();
            bestMatchDict.Add(bestMatch.Key, bestMatch.Value);

            return bestMatchDict;
        }

        /// <summary>
        /// returns non repeated matches and add to mods to ignore db combinations of mods that did not improve the matching.
        /// </summary>
        /// <param name="matchedPeptides"></param>
        /// <returns></returns>
        private Dictionary<PeptideWithSetModifications, Tuple<List<MatchedFragmentIon>, List<MatchedFragmentIon>>> AddNonCandidatesToModsToIgnore(
            IEnumerable<KeyValuePair<PeptideWithSetModifications, Tuple<List<MatchedFragmentIon>, List<MatchedFragmentIon>>>> matchedPeptides)
        {
            var sorted = matchedPeptides
                .OrderBy(x => x.Value.Item1[0].NeutralTheoreticalProduct.FragmentNumber)
                .ThenBy(x => x.Value.Item2[0].NeutralTheoreticalProduct.FragmentNumber);

            var goodMatches =
                new Dictionary<PeptideWithSetModifications,
                    Tuple<List<MatchedFragmentIon>, List<MatchedFragmentIon>>>();

            foreach (var match in sorted)
            {
                if (match.Value.Item1[0].NeutralTheoreticalProduct.FragmentNumber > sorted.First().Value.Item1[0].NeutralTheoreticalProduct.FragmentNumber ||
                    match.Value.Item2[0].NeutralTheoreticalProduct.FragmentNumber > sorted.First().Value.Item2[0].NeutralTheoreticalProduct.FragmentNumber)
                {
                    ModsToIgnore.Add(match.Key.AllModsOneIsNterminus);
                }
                else if (match.Value.Item1.Count >= sorted.First().Value.Item1.Count || match.Value.Item2.Count >= sorted.First().Value.Item2.Count)
                {
                    if (!goodMatches.Keys.Contains(match.Key))
                        goodMatches.Add(match.Key, match.Value);
                }
            }

            return goodMatches;
        }
        /// <summary>
        /// Returns a list of possible combination of mods that fit the delta mass provided. It assumes the Combination of Mods mass array is sorted,
        /// since it performs a binary search to get a range of acceptable mods.
        /// </summary>
        /// <param name="deltaMass"></param>
        /// <returns></returns>
        public List<List<Modification>> GetCombinationsThatFitDelta(double deltaMass)
        {
            var tolerance = new PpmTolerance(10);

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
        /// Recursive method that builds the combination of mods. 
        /// </summary>
        /// <param name="sortedListOfModsToAdd"></param>
        /// <param name="listOfModCombinations"></param>
        /// <param name="maxNumberOfModsInGroup"></param>
        /// <param name="allModsFromOneToN"></param>
        public static void CombinationBuilder(List<Modification> sortedListOfModsToAdd,
            ref List<List<Modification>> listOfModCombinations, int maxNumberOfModsInGroup, bool allModsFromOneToN)
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
                CombinationBuilder(sortedListOfModsToAdd, ref listOfModCombinations, maxNumberOfModsInGroup - 1, allModsFromOneToN);
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

        /// <summary>
        /// Returns the Modifications IdWithMotif joined. For CombinationBuilder use.
        /// </summary>
        /// <param name="list"></param>
        /// <returns></returns>
        private static string ModListNameString(List<Modification> list)
        {
            return String.Join("", list.Select(n => n.IdWithMotif));
        }

        /// <summary>
        /// Reads the PSMTSV file into a FilteredPsmTSV object (custom class for developing, later won't be necessary)
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
    }
}
