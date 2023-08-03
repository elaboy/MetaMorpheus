using System;
using EngineLayer;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common;
using iText.StyledXmlParser.Jsoup.Select;
using MzLibUtil;
using Proteomics.Fragmentation;
using TaskLayer;
using ThermoFisher.CommonCore.Data;
using UsefulProteomicsDatabases;
using Easy.Common.Extensions;
using FlashLFQ;

namespace Test
{
    public sealed class TestModComboSearch
    {
        private const string _filteredPsm = @"D:\08-30-22_bottomup\example.psmtsv";

        [Test]
        public void TestMultiModSearch()
        {
            var possibilities = MultiModSearch.GetPossibleModCombinations();
            var mods =
                Loaders.LoadUnimod(
                        @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                    .ToList();

            int maxToAdd = 3;
            var fixedMods = new List<Modification>();
            fixedMods.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")));

            var psmList =
                MultiModSearch.ReadFilteredPsmTSVShort(_filteredPsm);

            var msDataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML").GetAllScansList().Where(x => x.MsnOrder == 2).ToArray();
            var dataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML");
            List<Dictionary<Tuple<int, PeptideWithSetModifications>, List<MatchedFragmentIon>>> allMatches = new();

            foreach(var psm in psmList)
            {
                var peptideAsProtein = MultiModSearch.GetPeptideAsProtein(psm, fixedMods);

                var matches = MultiModSearch.GetPeptideFragmentIonsMatches(psm, dataFile, possibilities, fixedMods);
                allMatches.Add(matches);
            }
        }
    }

    public class MultiModSearch
    {
        public static Dictionary<Tuple<int,PeptideWithSetModifications>, List<MatchedFragmentIon>> GetPeptideFragmentIonsMatches(FilteredPsmTSV psm, MsDataFile dataFile,
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
                var ptm = GetPeptideWithMods(psm,  fixedMods, mod.ToList());
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
        private static double GetDeltaMass(FilteredPsmTSV psm, MsDataFile dataFile, PeptideWithSetModifications peptide)
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

            var deltaMass = GetDeltaMass(psm, dataFile, peptide.First());

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

            var maxIndex = ~ Array.BinarySearch(massArray, tolerance.GetMaximumValue(deltaMass));
            var minIndex = ~ Array.BinarySearch(massArray, tolerance.GetMinimumValue(deltaMass));

            int[] rangeIndex = Enumerable.Range(minIndex, Math.Abs(maxIndex-minIndex)).ToArray();

            var rangeOfPossibleMods = rangeIndex.Select(x => possibleMods.ToList()[x].Value).ToList();

            return rangeOfPossibleMods;
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

        private static IEnumerable<PeptideWithSetModifications> GetPeptideWithMods(FilteredPsmTSV psm, List<Modification> fixedMods = null, List<Modification> variableMods =null)
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

        private double GetDeltaMass()
        {
            return 0;
        }

        /// <summary>
        /// Returns all combinations of mods in the database. Currently hard-coded to three mods.
        /// </summary>
        /// <returns></returns>
        public static IOrderedEnumerable<KeyValuePair<double, Modification[]>> GetPossibleModCombinations(FilteredPsmTSV psm,
            MsDataFile dataFile, IEnumerable<PeptideWithSetModifications> peptide)
        {
            var tolerance = new PpmTolerance(40);

            var deltaMass = GetDeltaMass(psm, dataFile, peptide.First());

            var commonModsFromToml = GetModsFromGptmdThing().ToList();
            var mods =
                Loaders.LoadUnimod(
                        @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                    .ToList();

            commonModsFromToml.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C"))); //adds this mod into the toml database (it's a fixed mod)

            //var modsToAddPredicate = PredicateBuilder.Create<List<Modification>>(mod => mod);

            var groupedModsByOriginalId = from one in commonModsFromToml
                                          group one by one.OriginalId
                into newGroup
                                          select newGroup;

            var modCombos = (from grouped1 in groupedModsByOriginalId
                             from grouped2 in groupedModsByOriginalId
                             //from grouped3 in groupedModsByOriginalId
                             select new KeyValuePair<double, Modification[]>(
                                 (grouped1.First().MonoisotopicMass.Value + grouped2.First().MonoisotopicMass.Value),// + grouped3.First().MonoisotopicMass.Value),
                                 new[] { grouped1.First(), grouped2.First()})).OrderBy(x => x.Key);

            return modCombos;
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
