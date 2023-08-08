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
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Linq.Expressions;
using System.Reflection;
using System.Reflection.Metadata;
using System.Text.Json;
using System.Text.Json.Serialization;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Serialization;
using Easy.Common;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using pepXML.Generated;
using Proteomics.AminoAcidPolymer;
using TaskLayer;
using UsefulProteomicsDatabases;
using NetSerializer;
using OxyPlot;
using ScottPlot.Plottable;

namespace Test
{
    public sealed class TestModComboSearch
    {
        private const string _filteredPsm = @"D:\08-30-22_bottomup\example.psmtsv";

        [Test]
        public void TestMultiModSearch()
        {
            var mods =
                Loaders.LoadUnimod(
                        @"C:\Users\elabo\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                    .ToList();

            int maxToAdd = 3;
            var fixedMods = new List<Modification>();
            fixedMods.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")));

            var psmList =
                MultiModSearch.ReadFilteredPsmTSVShort(_filteredPsm);

            var msDataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML").GetAllScansList().Where(x => x.MsnOrder == 2).ToArray();
            var dataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML");
            List<KeyValuePair<Tuple<int, PeptideWithSetModifications>, List<MatchedFragmentIon>>> allMatches = new();

            for(int i = 1; i < 8;  i++)
            {
                foreach (var psm in psmList)
                {
                    var peptideAsProtein = MultiModSearch.GetPeptideAsProtein(psm, fixedMods);

                    var possibilities = MultiModSearch.GetPossibleModCombinations(psm, peptideAsProtein, i);


                    var matches = MultiModSearch.GetPeptideFragmentIonsMatches(psm, dataFile, possibilities, fixedMods);
                    var top = matches.OrderBy(x => x.Value.Count);
                    var selected = top.Last();
                    allMatches.Add(selected);
                }
            }

            var orderedMatches = allMatches.OrderBy(x => x.Value.Count()).GroupBy(x => x.Key.Item2.Protein.BaseSequence).ToList()
                .Select(x => new XMLGroup(){GroupName = x.Key, mods = new KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>(x.First().Key.Item2,x.First().Value)});

            


        }
    }

    public class XMLGroup
    {
        public string GroupName { get; set; }
        public KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>> mods { get; set; }
        
        internal XMLGroup() {}
    }

    public class MultiModSearch
    {
        public static Dictionary<Tuple<int, PeptideWithSetModifications>, List<MatchedFragmentIon>> GetPeptideFragmentIonsMatches(FilteredPsmTSV psm, MsDataFile dataFile,
            IOrderedEnumerable<KeyValuePair<double, Modification[]>> possibleMods, List<Modification> fixedMods)
        {
            Dictionary<Tuple<int, PeptideWithSetModifications>, List<MatchedFragmentIon>> results = new();

            var spectrum = dataFile.GetOneBasedScan(int.Parse(psm.ScanNumber));
            //var neutralMass = spectrum.SelectedIonMZ.Value.ToMass(spectrum.SelectedIonChargeStateGuess.Value);

            var probableMods =
                GetCombinationsThatFitDelta(possibleMods, GetPeptideWithMods(psm, fixedMods), psm, dataFile);

            var products = new List<Product>();
            int id = 0;
            foreach (var mod in probableMods)
            {
                var ptm = GetPeptideWithMods(psm, fixedMods, mod.ToList());
                foreach (var variant in ptm)
                {
                    variant.Fragment(spectrum.DissociationType ?? DissociationType.HCD, FragmentationTerminus.Both, products);

                    var match = MetaMorpheusEngine.MatchFragmentIons(new Ms2ScanWithSpecificMass(spectrum,
                        double.Parse(psm.PrecursorMass),
                        spectrum.SelectedIonChargeStateGuess.Value, dataFile.FilePath,
                        new CommonParameters()), products, new CommonParameters());

                    results.Add(new Tuple<int, PeptideWithSetModifications>(id, variant), match);
                    id = id + 1;
                }


            }
            return results;
        }



        private static IEnumerable<PeptideWithSetModifications> AddPTMs(FilteredPsmTSV psm, Modification[] mods, List<Modification> fixedMods)
        {
            var peptides = GetPeptideWithMods(psm, fixedMods, mods.ToList());

            return peptides.Select(peptide => new PeptideWithSetModifications(peptide.Protein, new DigestionParams(),
                oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: peptide.BaseSequence.Length,
                CleavageSpecificity.Full, "", 0, peptide.AllModsOneIsNterminus, peptide.NumFixedMods));

            //return new PeptideWithSetModifications(peptide.Protein, new DigestionParams(),
            //    oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: peptide.BaseSequence.Length,
            //    CleavageSpecificity.Full, "", 0, peptide.AllModsOneIsNterminus, peptide.NumFixedMods);
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
        /// Returns the Combination of mods that are possible. Uses binary search on the IOrderedEnumerable of mods. Currently hard-coded to three mods.
        /// </summary>
        /// <param name="possibleMods"></param>
        /// <param name="peptide"></param>
        /// <param name="psm"></param>
        /// <param name="dataFile"></param>
        /// <returns></returns>
        private static List<Modification[]> GetCombinationsThatFitDelta(IOrderedEnumerable<KeyValuePair<double, Modification[]>> possibleMods,
            IEnumerable<PeptideWithSetModifications> peptide, FilteredPsmTSV psm, MsDataFile dataFile)
        {
            var tolerance = new PpmTolerance(40);

            var deltaMass = GetDeltaMass(psm, peptide.First());

            var massArray = possibleMods.Select(x => x.Key).ToArray();

            var modsArray = possibleMods.Select(x => x.Value);

            var mods = from mod in modsArray
                       from mod1 in mod
                       from mod2 in mod
                           //from mod3 in mod
                       where psm.BaseSeq.Contains(mod1.Target.ToString()) && psm.BaseSeq.Contains(mod2.Target.ToString())// && psm.BaseSeq.Contains(mod3.Target.ToString())
                       select mod;

            var temp1 = tolerance.GetMaximumValue(deltaMass);
            var temp2 = tolerance.GetMinimumValue(deltaMass);

            var maxIndex = ~Array.BinarySearch(massArray, tolerance.GetMaximumValue(deltaMass));
            var minIndex = ~Array.BinarySearch(massArray, tolerance.GetMinimumValue(deltaMass));

            int[] rangeIndex = Enumerable.Range(minIndex, Math.Abs(maxIndex - minIndex)).ToArray();

            if (rangeIndex.Length == 0)
                rangeIndex = new []{maxIndex};
            var rangeOfPossibleMods = rangeIndex.Select(x => possibleMods.ToList()[x].Value).ToList();

            return rangeOfPossibleMods;
        }
        private static List<Modification[]> GetModificationsWithinDelta(
            IEnumerable<PeptideWithSetModifications> peptide, FilteredPsmTSV psm, MsDataFile dataFile, int numberOfVariableMods)
        {
            var tolerance = new PpmTolerance(40);


            var deltaMass = GetDeltaMass(psm, peptide.First());

            var commonModsFromToml = GetModsFromGptmdThing().ToList();
            var mods =
                Loaders.LoadUnimod(
                        @"C:\Users\elabo\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                    .ToList();

            commonModsFromToml.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C"))); //adds this mod into the toml database (it's a fixed mod)

            var groupedModsByOriginalId = (from one in commonModsFromToml
                                           group one by one.OriginalId
                into newGroup
                                           select newGroup).ToList();

            List<Modification[]> listOfMods = new();

            return new List<Modification[]>();
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

        /// <summary>
        /// Returns the the peptide without any mods from the given psm. Used for theoretical no mods mass.
        /// </summary>
        /// <param name="psm"></param>
        /// <returns></returns>
        private static IEnumerable<PeptideWithSetModifications> GetPeptideWithNoMods(FilteredPsmTSV psm)
        {
            var peptideWithoutMods = new Protein(psm.BaseSeq, psm.ProteinAccession)
                .Digest(new DigestionParams("top-down"), new List<Modification>(), new List<Modification>());

            return peptideWithoutMods;
        }

        /// <summary>
        /// Returns the peptide from given psm with fixed modifications.
        /// </summary>
        /// <param name="digestedProtein"></param>
        /// <param name="psm"></param>
        /// <returns></returns>
        private static PeptideWithSetModifications GetPeptideWithFixedMods(IEnumerable<PeptideWithSetModifications> digestedProtein, FilteredPsmTSV psm)
        {
            var unmodifiedPeptide = new PeptideWithSetModifications(digestedProtein.First().Protein,
                new DigestionParams(), 1, digestedProtein.First().BaseSequence.Length,
                CleavageSpecificity.Full, "",
                0, digestedProtein.First().AllModsOneIsNterminus,
                digestedProtein.First().NumFixedMods);

            return unmodifiedPeptide;
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

        public static IEnumerable<KeyValuePair<double, Modification[]>> GetKCombinations(
            IEnumerable<IGrouping<string, Modification>> mods, int numberOfMods)
        {
            var combination = new List<KeyValuePair<double, Modification[]>>();

            switch (numberOfMods)
            {
                case 1:
                    var oneMod = from mod in mods
                        select new KeyValuePair<double, Modification[]>(
                            mod.First().MonoisotopicMass.Value, new[] { mod.First() });
                    combination = oneMod.ToList();
                    break;

                case 2:
                    var twoMods = from mod in mods
                        from mod2 in mods
                        select new KeyValuePair<double, Modification[]>(
                            mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value,
                            new[] { mod.First(), mod2.First() });
                    combination = twoMods.ToList();
                    break;

                case 3:
                    var threeMods = from mod in mods
                        from mod2 in mods
                        from mod3 in mods
                        select new KeyValuePair<double, Modification[]>(
                            mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value + mod3.First().MonoisotopicMass.Value,
                            new[] { mod.First(), mod2.First(), mod3.First() });
                    combination = threeMods.ToList();
                    break;

                case 4:
                    var fourMods = from mod in mods
                        from mod2 in mods
                        from mod3 in mods
                        from mod4 in mods
                        select new KeyValuePair<double, Modification[]>(
                            mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value + mod3.First().MonoisotopicMass.Value
                            + mod4.First().MonoisotopicMass.Value,
                            new[] { mod.First(), mod2.First(), mod3.First(), mod4.First() });
                    combination = fourMods.ToList();
                    break;

                case 5:
                    var fiveMods = from mod in mods
                        from mod2 in mods
                        from mod3 in mods
                        from mod4 in mods
                        from mod5 in mods
                        select new KeyValuePair<double, Modification[]>(
                            mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value + mod3.First().MonoisotopicMass.Value
                            + mod4.First().MonoisotopicMass.Value + mod5.First().MonoisotopicMass.Value,
                            new[] { mod.First(), mod2.First(), mod3.First(), mod4.First(), mod5.First() });
                    combination = fiveMods.ToList();
                    break;

                case 6:
                    var sixMods = from mod in mods
                        from mod2 in mods
                        from mod3 in mods
                        from mod4 in mods
                        from mod5 in mods
                        from mod6 in mods
                        select new KeyValuePair<double, Modification[]>(
                            mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value + mod3.First().MonoisotopicMass.Value
                            + mod4.First().MonoisotopicMass.Value + mod5.First().MonoisotopicMass.Value + mod6.First().MonoisotopicMass.Value,
                            new[] { mod.First(), mod2.First(), mod3.First(), mod4.First(), mod5.First(), mod6.First() });
                    combination = sixMods.ToList();
                    break;

                case 7:
                    var sevenMods = from mod in mods
                        from mod2 in mods
                        from mod3 in mods
                        from mod4 in mods
                        from mod5 in mods
                        from mod6 in mods
                        from mod7 in mods
                        select new KeyValuePair<double, Modification[]>(
                            mod.First().MonoisotopicMass.Value + mod2.First().MonoisotopicMass.Value + mod3.First().MonoisotopicMass.Value
                            + mod4.First().MonoisotopicMass.Value + mod5.First().MonoisotopicMass.Value + mod6.First().MonoisotopicMass.Value
                            + mod7.First().MonoisotopicMass.Value,
                            new[] { mod.First(), mod2.First(), mod3.First(), mod4.First(), mod5.First(), mod6.First(), mod7.First() });
                    combination = sevenMods.ToList();
                    break;
            }

            return combination.AsEnumerable();
        }
        /// <summary>
        /// Returns all combinations of mods in the database. Currently hard-coded to three mods.
        /// </summary>
        /// <returns></returns>
        public static IOrderedEnumerable<KeyValuePair<double, Modification[]>> GetPossibleModCombinations(FilteredPsmTSV psm,
                IEnumerable<PeptideWithSetModifications> peptide, int numberOfVariableMods)
        {
            //var tolerance = new PpmTolerance(40);

            //var deltaMass = GetDeltaMass(psm, dataFile, peptide.First());

            var commonModsFromToml = GetModsFromGptmdThing()
                .Where(mod => mod.Target.ToString().Contains('X') ||
                              peptide.First().BaseSequence.Contains(mod.Target.ToString())).ToList();

            var mods =
                Loaders.LoadUnimod(
                        @"C:\Users\elabo\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                    .ToList();

            commonModsFromToml.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C"))); //adds this mod into the toml database (it's a fixed mod)

            var deltaMass = GetDeltaMass(psm, peptide.First());

            var modsGroupedByOriginalId = commonModsFromToml.GroupBy(modId => modId.IdWithMotif);

            return GetKCombinations(modsGroupedByOriginalId, 3).OrderBy(x => x.Key);
        }

        /// <summary>
        /// Returns mods from a toml containing common modifications. 
        /// </summary>
        /// <param name="gptmdToml"></param>
        /// <returns></returns>
        private static IEnumerable<Modification> GetModsFromGptmdThing(string gptmdToml = @"Task1-GPTMDTaskconfig.toml")
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
    }
}
