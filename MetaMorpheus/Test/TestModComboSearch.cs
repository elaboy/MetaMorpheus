using Easy.Common.Extensions;
using EngineLayer;
using EngineLayer.CombinatorialSearch;
using ExCSS;
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
using System.Data;
using System.IO;
using System.Linq;
using System.Text.Json;
using System.Threading.Tasks;
using TaskLayer;
using TaskLayer.CombinatorialSearchTask;
using ThermoFisher.CommonCore.Data;
using UsefulProteomicsDatabases;

namespace Test
{
    public sealed class TestModComboSearch
    {
        private const string _filteredPsm = @"D:\08-30-22_bottomup\example.psmtsv";

        #region TEST ENCODING MANUAL

        //MsDataFile
        private const string _dataFile =
            @"D:\08-30-22_bottomup\fractionated\08-31-22_fractionated_human_Tryp_40ug_F7.raw";

        //peptide 1
        private string[] _peptide1 = new[]
        {
            "08-31-22_fractionated_human_Tryp_40ug_F7-calib",
            "25726",
            "25711",
            "44.376",
            "AASAAAASAAAASAASGSPGPGEGSAGGEKR",
            "[UniProt:N-acetylalanine on A]AASAAAASAAAASAAS[Common Biological:Phosphorylation on S]GSPGPGEGSAGGEKR",
            "N-acetylalanine on A Phosphorylation on S",
            "2",
            "Q13263",
            "Transcription intermediary factor 1-beta",
            "primary:TRIM28, synonym:KAP1, synonym:RNF96, synonym:TIF1B",
            "Homo sapiens",
            "2 to 32",
            "[y1+1, y2+1, y4+1, y5+1, y6+1, y7+1, y8+1, y9+1, y10+1, y11+1, y12+1, y13+1, y16+1, (y16-97.98)+1, y17+1, (y17-97.98)+1, y18+1, (y18-97.98)+1, y19+1, (y19-97.98)+1, y20+1, (y20-97.98)+1, y21+1, (y21-97.98)+1, y22+1, (y22-97.98)+1, (y23-97.98)+1];[b2+1, b3+1, b4+1, b5+1, b6+1, b7+1, b8+1, b9+1, b10+1, b11+1, b12+1, b13+1, b15+1, (b24-97.98)+1];[(M0-97.98)+2]",
            "42",
            "2664.17846",
            "2"
        };

        //peptide 2

        private string[] _peptide2 = new[]
        {
            "08-31-22_fractionated_human_Tryp_40ug_F7-calib",
            "49967",
            "49954",
            "38.301",
            "LASPSGSTSSGLEVVAPEGTSAPGGGPGTLDDSATICR",
            "LASPS[Common Biological:Phosphorylation on S]GSTSSGLEVVAPEGTSAPGGGPGTLDDSATIC[Common Fixed:Carbamidomethyl on C]R",
            "Carbamidomethyl on C Phosphorylation on S",
            "1",
            "Q13263",
            "Transcription intermediary factor 1-beta",
            "primary:TRIM28, synonym:KAP1, synonym:RNF96, synonym:TIF1B",
            "Homo sapiens",
            "592 to 629",
            "[y1+1, y2+1, y3+1, y4+1, y5+1, y6+1, y7+1, y8+1, y9+1, y10+1, y11+1, y12+1, y13+1, y14+1, y15+1, y16+2, y17+1, y18+1, y19+1, y20+1, y22+2, y23+2];[(b6-97.98)+1, (b7-97.98)+1, (b8-97.98)+1, (b10-97.98)+1, b11+1, (b11-97.98)+1, b12+1, (b12-97.98)+1, b13+1, (b13-97.98)+1, b14+1, (b14-97.98)+1, b15+1, (b15-97.98)+1, b16+1, (b16-97.98)+1]",
            "38",
            "3637.64432",
            "3"
        };

        #endregion

        [Test]
        public void TestEncodingMatchingManuallyTwoExamples()
        {
            var mods =
                Loaders.LoadUnimod(
                        @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                    .ToList();

            var fixedMods = new List<Modification>();
            fixedMods.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")));


            var dataFile = MsDataFileReader.GetDataFile(_dataFile);
            var psm1 = new FilteredPsmTSV()
            {
                FileName = _peptide1[0],
                ScanNumber = _peptide1[1],
                PrecursorScanNumber = _peptide1[2],
                Score = _peptide1[3],
                BaseSeq = _peptide1[4],
                FullSeq = _peptide1[5],
                Mods = _peptide1[6],
                ModsCount = _peptide1[7],
                ProteinAccession = _peptide1[8],
                ProteinName = _peptide1[9],
                GeneName = _peptide1[10],
                OrganismName = _peptide1[11],
                StartAndEndResiduesInProtein = _peptide1[12],
                MatchedIonSeries = _peptide1[13],
                MatchedIonCounts = _peptide1[14],
                PrecursorMass = _peptide1[15],
                Charge = _peptide1[16]
            };
            var scanForPsm1 = dataFile.GetOneBasedScan(int.Parse(psm1.ScanNumber));
            var psm2 = new FilteredPsmTSV(_peptide2);

            //List<Modification> commonBiologycalMods = GlobalVariables.AllModsKnown.OfType<Modification>()
            //    .Where(mod => mod.IdWithMotif.Contains("Phosphorylation on S") || mod.IdWithMotif.Contains("N-acetylalanine")).ToList();

            List<Modification> commonBiologycalMods = GlobalVariables.AllModsKnown.OfType<Modification>()
                .Where(mod => mod.ModificationType.Contains("Common Biological") ||
                              mod.ModificationType.Contains("Common Artifact") ||
                              mod.ModificationType.Contains("Less Common") ||
                              mod.ModificationType.Contains("Metals") ||
                              mod.ModificationType.Contains("UniProt")).ToList();

            //List<Modification> commonBiologycalMods = GlobalVariables.AllModsKnown.OfType<Modification>()
            //    .Where(mod => mod.IdWithMotif.Contains("Phosphorylation on S") || mod.IdWithMotif.Contains("N-acetylalanine")).ToList();

            var engine = new MultipleSearchEngine(new List<FilteredPsmTSV>() { psm1 }, commonBiologycalMods, 3,
                fixedMods, dataFile, true);

            List<Tuple<PeptideWithSetModifications, List<MatchedFragmentIon>>> tempMatches = new();

            var protein = new Protein(psm1.BaseSeq, psm1.ProteinAccession);

            var peptideFromProtein = new PeptideWithSetModifications(protein, new DigestionParams("top-down"), 1,
                psm1.BaseSeq.Length,
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 1);

            var fixedModPeptide = peptideFromProtein.Protein.Digest(new DigestionParams("top-down", 0), fixedMods,
                new List<Modification>());

            var deltaMass = Math.Abs(fixedModPeptide.First().MonoisotopicMass - double.Parse(psm1.PrecursorMass));

            var possibleComboMods = engine.GetCombinationsThatFitDelta(deltaMass);

            var bProducts = new List<Product>();
            var yProducts = new List<Product>();

            Dictionary<PeptideWithSetModifications, Tuple<List<MatchedFragmentIon>, List<MatchedFragmentIon>>>
                matchesList = new();

            List<Dictionary<int, Modification>> modsToIgnore = new(); //mods to ignore, to avoid unnecessary checks 


            foreach (var combo in possibleComboMods)
            {
                var comboToTry =
                    peptideFromProtein.Protein.Digest(new DigestionParams("top-down", 0), fixedMods, combo);

                Dictionary<PeptideWithSetModifications, Tuple<List<MatchedFragmentIon>, List<MatchedFragmentIon>>>
                    tempMatchesFragmentIons = new(); //b and y fragments tuple

                List<PeptideWithSetModifications> filteredPeptides = new();
                if (modsToIgnore.Count > 0)
                {
                    foreach (var peptide in comboToTry)
                    {
                        if (!modsToIgnore.Contains(peptide.AllModsOneIsNterminus))
                        {
                            filteredPeptides.Add(peptide);
                        }
                    }
                }
                else
                {
                    {
                        filteredPeptides = comboToTry.ToList();
                    }
                }


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
                        tempMatchesFragmentIons.Add(peptide,
                            new Tuple<List<MatchedFragmentIon>, List<MatchedFragmentIon>>(bMatch, yMatch));

                }

                var tempMatchesFragmentIonsWithoutZeros =
                    tempMatchesFragmentIons.Where(x => x.Value.Item1.Count > 0 && x.Value.Item2.Count > 0);

                var sorted = tempMatchesFragmentIonsWithoutZeros
                    .OrderBy(x => x.Value.Item1[0].NeutralTheoreticalProduct.FragmentNumber)
                    .ThenBy(x => x.Value.Item2[0].NeutralTheoreticalProduct.FragmentNumber);

                //Get rid of bad variants
                foreach (var match in sorted)
                {
                    if (match.Value.Item1[0].NeutralTheoreticalProduct.FragmentNumber >
                        sorted.First().Value.Item1[0].NeutralTheoreticalProduct.FragmentNumber ||
                        match.Value.Item2[0].NeutralTheoreticalProduct.FragmentNumber > sorted.First().Value.Item2[0]
                            .NeutralTheoreticalProduct.FragmentNumber)
                    {
                        modsToIgnore.Add(match.Key.AllModsOneIsNterminus);
                    }
                    else if (match.Value.Item1.Count >= sorted.First().Value.Item1.Count ||
                             match.Value.Item2.Count >= sorted.First().Value.Item2.Count)
                    {
                        if (!matchesList.Keys.Contains(match.Key))
                            matchesList.Add(match.Key, match.Value);
                    }
                }
            }

            var resultsBByDescending = matchesList.OrderByDescending(x => x.Value.Item1.Count);
            var resultsYByDescending = matchesList.OrderByDescending(x => x.Value.Item2.Count);
            var groupedNoGo = modsToIgnore.SelectMany(x => x.Values);
            var bestModsFromResultsBandY = new Dictionary<int, Modification>();
            bestModsFromResultsBandY.Add(resultsBByDescending.First().Key.AllModsOneIsNterminus);
            bestModsFromResultsBandY.Add(resultsYByDescending.First().Key.AllModsOneIsNterminus);


            var secondRunPeptide = new PeptideWithSetModifications(resultsYByDescending.First().Key.Protein,
                new DigestionParams("top-down", 0), 1, resultsYByDescending.First().Key.BaseSequence.Length,
                CleavageSpecificity.Full, "", 0, bestModsFromResultsBandY, 1);

            var newDeltaMass = secondRunPeptide.MonoisotopicMass -
                               double.Parse(psm1.PrecursorMass);

            var secondRunDigested = secondRunPeptide.Protein.Digest(new DigestionParams("top-down", 0), fixedMods,
                bestModsFromResultsBandY.Values.ToList());

            var secondProductsRun = new List<Product>();
            var secondRunResults = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();
            foreach (var product in secondRunDigested)
            {
                product.Fragment(DissociationType.HCD, FragmentationTerminus.Both, secondProductsRun);

                var secondRunMatch = MetaMorpheusEngine.MatchFragmentIons(
                    new Ms2ScanWithSpecificMass(scanForPsm1, scanForPsm1.SelectedIonMZ.Value,
                        int.Parse(psm1.Charge), dataFile.FilePath,
                        new CommonParameters()), secondProductsRun, new CommonParameters());

                secondRunResults.Add(product, secondRunMatch);
            }

            var finalMatchResult = secondRunResults.OrderByDescending(x => x.Value.Count).First();
            Assert.Equals(true, false);

        }

        [Test]
        public void TestEncodingMatchingManuallyFixedModAsVariableTry()
        {
            var mods =
                Loaders.LoadUnimod(
                        @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                    .ToList();

            var fixedMods = new List<Modification>();


            var dataFile = MsDataFileReader.GetDataFile(_dataFile);
            var psm1 = new FilteredPsmTSV()
            {
                FileName = _peptide1[0],
                ScanNumber = _peptide1[1],
                PrecursorScanNumber = _peptide1[2],
                Score = _peptide1[3],
                BaseSeq = _peptide1[4],
                FullSeq = _peptide1[5],
                Mods = _peptide1[6],
                ModsCount = _peptide1[7],
                ProteinAccession = _peptide1[8],
                ProteinName = _peptide1[9],
                GeneName = _peptide1[10],
                OrganismName = _peptide1[11],
                StartAndEndResiduesInProtein = _peptide1[12],
                MatchedIonSeries = _peptide1[13],
                MatchedIonCounts = _peptide1[14],
                PrecursorMass = _peptide1[15],
                Charge = _peptide1[16]
            };
            var scanForPsm1 = dataFile.GetOneBasedScan(int.Parse(psm1.ScanNumber));
            var psm2 = new FilteredPsmTSV(_peptide2);

            //List<Modification> commonBiologycalMods = GlobalVariables.AllModsKnown.OfType<Modification>()
            //    .Where(mod => mod.IdWithMotif.Contains("Phosphorylation on S") || mod.IdWithMotif.Contains("N-acetylalanine")).ToList();

            List<Modification> commonBiologycalMods = GlobalVariables.AllModsKnown.OfType<Modification>()
                .Where(mod => mod.ModificationType.Contains("Common Biological") ||
                              mod.ModificationType.Contains("Common Artifact") ||
                              mod.ModificationType.Contains("Less Common") ||
                              mod.ModificationType.Contains("Metals") ||
                              mod.ModificationType.Contains("UniProt")).ToList();

            commonBiologycalMods.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")));


            //List<Modification> commonBiologycalMods = GlobalVariables.AllModsKnown.OfType<Modification>()
            //    .Where(mod => mod.IdWithMotif.Contains("Phosphorylation on S") || mod.IdWithMotif.Contains("N-acetylalanine")).ToList();

            var engine = new MultipleSearchEngine(new List<FilteredPsmTSV>() { psm1 }, commonBiologycalMods, 3,
                fixedMods, dataFile, true);

            List<Tuple<PeptideWithSetModifications, List<MatchedFragmentIon>>> tempMatches = new();

            var protein = new Protein(psm1.BaseSeq, psm1.ProteinAccession);

            var peptideFromProtein = new PeptideWithSetModifications(protein, new DigestionParams("top-down"), 1,
                psm1.BaseSeq.Length,
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 1);

            var fixedModPeptide = peptideFromProtein.Protein.Digest(new DigestionParams("top-down", 0), fixedMods,
                new List<Modification>());

            var deltaMass = Math.Abs(fixedModPeptide.First().MonoisotopicMass - double.Parse(psm1.PrecursorMass));

            var possibleComboMods = engine.GetCombinationsThatFitDelta(deltaMass);

            var bProducts = new List<Product>();
            var yProducts = new List<Product>();

            Dictionary<PeptideWithSetModifications, Tuple<List<MatchedFragmentIon>, List<MatchedFragmentIon>>>
                matchesList = new();

            List<Dictionary<int, Modification>> modsToIgnore = new(); //mods to ignore, to avoid unnecessary checks 


            foreach (var combo in possibleComboMods)
            {
                var comboToTry =
                    peptideFromProtein.Protein.Digest(new DigestionParams("top-down", 0), fixedMods, combo);

                Dictionary<PeptideWithSetModifications, Tuple<List<MatchedFragmentIon>, List<MatchedFragmentIon>>>
                    tempMatchesFragmentIons = new(); //b and y fragments tuple

                List<PeptideWithSetModifications> filteredPeptides = new();
                if (modsToIgnore.Count > 0)
                {
                    foreach (var peptide in comboToTry)
                    {
                        if (!modsToIgnore.Contains(peptide.AllModsOneIsNterminus))
                        {
                            filteredPeptides.Add(peptide);
                        }
                    }
                }
                else
                {
                    {
                        filteredPeptides = comboToTry.ToList();
                    }
                }


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
                        tempMatchesFragmentIons.Add(peptide,
                            new Tuple<List<MatchedFragmentIon>, List<MatchedFragmentIon>>(bMatch, yMatch));

                }

                var tempMatchesFragmentIonsWithoutZeros =
                    tempMatchesFragmentIons.Where(x => x.Value.Item1.Count > 0 && x.Value.Item2.Count > 0);

                var sorted = tempMatchesFragmentIonsWithoutZeros
                    .OrderBy(x => x.Value.Item1[0].NeutralTheoreticalProduct.FragmentNumber)
                    .ThenBy(x => x.Value.Item2[0].NeutralTheoreticalProduct.FragmentNumber);

                //Get rid of bad variants
                foreach (var match in sorted)
                {
                    if (match.Value.Item1[0].NeutralTheoreticalProduct.FragmentNumber >
                        sorted.First().Value.Item1[0].NeutralTheoreticalProduct.FragmentNumber ||
                        match.Value.Item2[0].NeutralTheoreticalProduct.FragmentNumber > sorted.First().Value.Item2[0]
                            .NeutralTheoreticalProduct.FragmentNumber)
                    {
                        modsToIgnore.Add(match.Key.AllModsOneIsNterminus);
                    }
                    else if (match.Value.Item1.Count >= sorted.First().Value.Item1.Count ||
                             match.Value.Item2.Count >= sorted.First().Value.Item2.Count)
                    {
                        if (!matchesList.Keys.Contains(match.Key))
                            matchesList.Add(match.Key, match.Value);
                    }
                }
            }

            var resultsBByDescending = matchesList.OrderByDescending(x => x.Value.Item1.Count);
            var resultsYByDescending = matchesList.OrderByDescending(x => x.Value.Item2.Count);
            var groupedNoGo = modsToIgnore.SelectMany(x => x.Values);
            var bestModsFromResultsBandY = new Dictionary<int, Modification>();
            bestModsFromResultsBandY.Add(resultsBByDescending.First().Key.AllModsOneIsNterminus);
            bestModsFromResultsBandY.Add(resultsYByDescending.First().Key.AllModsOneIsNterminus);


            var secondRunPeptide = new PeptideWithSetModifications(resultsYByDescending.First().Key.Protein,
                new DigestionParams("top-down", 0), 1, resultsYByDescending.First().Key.BaseSequence.Length,
                CleavageSpecificity.Full, "", 0, bestModsFromResultsBandY, 1);

            var newDeltaMass = secondRunPeptide.MonoisotopicMass -
                               double.Parse(psm1.PrecursorMass);

            var secondRunDigested = secondRunPeptide.Protein.Digest(new DigestionParams("top-down", 0), fixedMods,
                bestModsFromResultsBandY.Values.ToList());

            var secondProductsRun = new List<Product>();
            var secondRunResults = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();
            foreach (var product in secondRunDigested)
            {
                product.Fragment(DissociationType.HCD, FragmentationTerminus.Both, secondProductsRun);

                var secondRunMatch = MetaMorpheusEngine.MatchFragmentIons(
                    new Ms2ScanWithSpecificMass(scanForPsm1, scanForPsm1.SelectedIonMZ.Value,
                        int.Parse(psm1.Charge), dataFile.FilePath,
                        new CommonParameters()), secondProductsRun, new CommonParameters());

                secondRunResults.Add(product, secondRunMatch);
            }

            var finalMatchResult = secondRunResults.OrderByDescending(x => x.Value.Count).First();
            Assert.Equals(true, false);

        }

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
            File.WriteAllLines(@"C:\Users\Edwin\Desktop\mods.txt",
                completeListOfModCombinations.Select(n => ModListNameString(n)).ToArray());
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


            List<DataTable> proteinGroupsTables = new();

            Parallel.ForEach(results, result =>
            {
                var table = new DataTable();
                foreach (var feature in typeof(MultiModSearchResults).GetProperties())
                {
                    table.Columns.Add(new DataColumn(feature.Name));
                }

                foreach (var peptide in result)
                {
                    var individualPeptide = new MultiModSearchResults()
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
                    };
                    DataRow row = table.NewRow();
                    row[1] = individualPeptide.BaseSequece;
                    row[2] = individualPeptide.AccessionNumber;
                    row[3] = individualPeptide.PeptideLength;
                    row[3] = individualPeptide.MonoisotopicMass;
                    row[4] = individualPeptide.MostAbundantMonoisotopicMass;
                    row[5] = individualPeptide.IsDecoy;
                    row[6] = String.Join(", ", individualPeptide.MatchedIons);
                    row[7] = String.Join(", ", individualPeptide.MatchedIonCharge);
                    row[8] = String.Join(", ", individualPeptide.TheoricalMz);
                    row[9] = String.Join(", ", individualPeptide.MatchedMz);
                    row[10] = String.Join(", ", individualPeptide.MassErrorPpm);
                    row[11] = String.Join(", ", individualPeptide.MassErrorDa);
                }

                proteinGroupsTables.Add(table);
            });


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
        public void METHOD()
        {
            var temp = new List<PsmFromTsv>();


            temp.Add(PsmTsvReader.ReadTsv(
                @"D:\topDown\MOxAndBioMetArtModsGPTMD_Search\Task2-SearchTask\AllPSMs.psmtsv",
                out List<string> warnings));
            var baseSeqGroup = temp.GroupBy(x => new { x.BaseSeq });
            var temp2 = temp.GroupBy(p => new
            {
                p.FullSequence,
                p.PreviousAminoAcid,
                p.NextAminoAcid
            });
        }

        [Test]
        public static void Bubba()
        {
            List<Modification> allAvailableMods = PtmListLoader
                .ReadModsFromFile(
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "CommonBiological.txt"),
                    out var errors).ToList();

            List<Modification> chosenSubsetOfMods = allAvailableMods.ToArray()[0..5].ToList();
            List<List<Modification>> completeListOfModCombinations = new();
            Mc(chosenSubsetOfMods, ref completeListOfModCombinations, 3, true);
            File.WriteAllLines(@"E:\junk\mods.txt",
                completeListOfModCombinations.Select(n => ModListNameString(n)).ToArray());
            Assert.AreEqual(1, 0);
        }

        public static void Mc(List<Modification> sortedListOfModsToAdd,
            ref List<List<Modification>> listOfModCombinations,
            int maxNumberOfModsInGroup, bool allModsFromOneToN)
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
                    listOfModCombinations =
                        listOfModCombinations.Where(c => c.Count == maxNumberOfModsInGroup).ToList();
                }
            }

        }

        public static string ModListNameString(List<Modification> list)
        {
            return String.Join("", list.Select(n => n.IdWithMotif));
        }

        [Test]
        public void METHOD23()
        {
            var taskList = new List<(string, MetaMorpheusTask)>();
            var CSTask = new CSTask(new CommonParameters());

            taskList.Add(("CS-Task", CSTask));

            List<DbForTask> dbForTask = new List<DbForTask>();
            dbForTask.Add(new DbForTask(
                @"D:\08-30-22_bottomup\uniprotkb_accession_A0A0C5B5G6_OR_access_2023_09_01.fasta",
                false));

            var runner = new EverythingRunnerEngine(taskList,
                new List<string>() { @"D:\08-30-22_bottomup\test.mzML" },
                dbForTask, @"D:\TestingCSTask");

            runner.Run();
        }

        [Test]
        public void METHOD23_FullFiles()
        {
            var taskList = new List<(string, MetaMorpheusTask)>();
            var CSTask = new CSTask(new CommonParameters());

            taskList.Add(("CS-Task", CSTask));

            List<DbForTask> dbForTask = new List<DbForTask>();
            dbForTask.Add(new DbForTask(
                @"D:\08-30-22_bottomup\uniprotkb_taxonomy_id_9606_AND_reviewed_2023_09_04.fasta",
                false));

            var runner = new EverythingRunnerEngine(taskList,
                new List<string>()
                {
                    @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F8-calib.mzML",
                    @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F1-calib.mzML",
                    @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F2-calib.mzML",
                    @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F3-calib.mzML",
                    @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F4-calib.mzML",
                    @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F5-calib.mzML",
                    @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F6-calib.mzML",
                    @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F7-calib.mzML"
                },
                dbForTask, @"D:\TestingCSTask");

            runner.Run();
        }

        [Test]
        public void METHOD23_Unfractionated()
        {
            var taskList = new List<(string, MetaMorpheusTask)>();
            var CSTask = new CSTask(new CommonParameters());

            taskList.Add(("CS-Task", CSTask));

            List<DbForTask> dbForTask = new List<DbForTask>();
            dbForTask.Add(new DbForTask(
                @"D:\08-30-22_bottomup\uniprotkb_taxonomy_id_9606_AND_reviewed_2023_09_04.fasta",
                false));

            var runner = new EverythingRunnerEngine(taskList,
                new List<string>()
                {
                    @"D:\08-30-22_bottomup\UnfractionatedCalibrated\Task1-CalibrateTask\08-29-22_unfractionated_human_Tryp-calib.mzML"
                },
                dbForTask, @"D:\TestingCSTask");

            runner.Run();

            var xmlProteins = UsefulProteomicsDatabases.ProteinDbLoader.LoadProteinXML(
                @"D:\08-29-22_unfractionated_human_Tryp_TestingCSSearch\08-29-22_unfractionated_human_Tryp.xml", true, DecoyType.None,
                GlobalVariables.AllModsKnown.OfType<Modification>()
                    .Where(mod =>
                        mod.ModificationType.Contains("Common Biological") ||
                        mod.ModificationType.Contains("Common Artifact") ||
                        mod.ModificationType.Contains("Less Common") ||
                        mod.ModificationType.Contains("Metals") ||
                        mod.ModificationType.Contains("UniProt")), false, new[] { "null" }, out _);

            var proteinsWithMods = xmlProteins.Where(x => x.OneBasedPossibleLocalizedModifications.Count > 0).ToList()
                .OrderByDescending(x => x.OneBasedPossibleLocalizedModifications.Count);

        }

        [Test]
        public void LOAD_XML()
        {
            var xmlProteins = UsefulProteomicsDatabases.ProteinDbLoader.LoadProteinXML(
                @"D:\TestingCSTask\testingTasks.xml", true, DecoyType.None,
                GlobalVariables.AllModsKnown.OfType<Modification>()
                    .Where(mod =>
                        mod.ModificationType.Contains("Common Biological") ||
                        mod.ModificationType.Contains("Common Artifact") ||
                        mod.ModificationType.Contains("Less Common") ||
                        mod.ModificationType.Contains("Metals") ||
                        mod.ModificationType.Contains("UniProt")), false, new[] { "null" }, out _);

            var proteinsWithMods = xmlProteins.Where(x => x.OneBasedPossibleLocalizedModifications.Count > 0).ToList()
                .OrderByDescending(x => x.OneBasedPossibleLocalizedModifications.Count);
        }

        [Test]
        public void GetAllAccessionNumbersAndFeaturesCSTask()
        {
            var taskList = new List<(string, MetaMorpheusTask)>();
            var CSTask = new CSTask(new CommonParameters());

            taskList.Add(("CS-Task", CSTask));

            List<DbForTask> dbForTask = new List<DbForTask>();
            dbForTask.Add(new DbForTask(
                @"D:\08-30-22_bottomup\uniprotkb_proteome_UP000005640_AND_revi_2023_09_05.xml",
                false));

            var runner = new EverythingRunnerEngine(taskList,
                new List<string>() {
                    @"D:\08-30-22_bottomup\test.mzML" },
                dbForTask, @"D:\TestingCSTask");

            runner.Run();
        }

        [Test]
        public void TestRefactoresSearchEngineOnTopDownData()
        {
            var mods =
                Loaders.LoadUnimod(
                        @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml").ToList();

            var fixedMods = new List<Modification>();
            //fixedMods.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C"))); toml indicates no fixed mods


            var dataFile = MsDataFileReader.GetDataFile(@"D:\topDown\test.mzML");

            List<Modification> commonBiologycalMods = GlobalVariables.AllModsKnown.OfType<Modification>()
                .Where(mod =>
                    mod.ModificationType.Contains("Common Biological") ||
                    mod.ModificationType.Contains("Common Artifact") ||
                    mod.ModificationType.Contains("Less Common") ||
                    mod.ModificationType.Contains("Metals") ||
                    mod.ModificationType.Contains("UniProt")).ToList();

            var psms = MMGPTMD.ReadFilteredPsmTSVShort(@"D:\topDown\example.psmtsv");

            var msDataFile = Readers.MsDataFileReader.GetDataFile(@"D:\topDown\test.mzML").LoadAllStaticData();

            var engine = new MultipleSearchEngine(psms, commonBiologycalMods, 3,
                fixedMods, msDataFile, true);

            var results = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();

            foreach (var protein in engine.ProteinListInferedFromGPTMD)
            {
                var peptide = protein.Item2
                    .Digest(new DigestionParams("top-down", 0), fixedMods,
                    new List<Modification>());

                var deltaMass = protein.Item1 - peptide.First().MonoisotopicMass;

                var possibleMods = engine.GetCombinationsThatFitDelta(deltaMass);

                var psm = psms.Find(x => x.BaseSeq.Equals(protein.Item2.BaseSequence));

                foreach (var mod in possibleMods)
                {
                    var bestCandidate = engine.SearchForMods(fixedMods, possibleMods,
                        protein.Item2.Digest(new DigestionParams("top-down", 0),
                            fixedMods, mod),
                        msDataFile.GetOneBasedScan(int.Parse(psm.ScanNumber)), psm, msDataFile);

                    if (results.Count == 0)
                    {
                        results.Add(bestCandidate);
                        continue;
                    }

                    if (bestCandidate.Count == 0)
                    {
                        continue;
                    }

                    if (!results.ContainsKey(bestCandidate.Keys.First()))
                    {
                        results.Add(bestCandidate);
                    }
                }


            }

            engine.WriteToXMLDatabase(@"D:\topDown\databaseFromMods.xml",
                engine.ProteinListInferedFromGPTMD.Select(x => x.Item2).ToList(),
                new List<string>() { msDataFile.FilePath }, results);
            //var groupedResults = results.GroupBy(x => x.Key.BaseSequence);
        }
    }


    public class Comparing
    {
        public string peptideSequence { get; set; }
        public int mmMatchCount { get; set; }
        public int myMatchCount { get; set; }

    }

    public class MultiModSearch
    {
        public KeyValuePair<double, Modification[]>[] CombinationsFromDatabase { get; set; }
        public IOrderedEnumerable<IGrouping<string, Modification>> ModificationGroups { get; set; }
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
