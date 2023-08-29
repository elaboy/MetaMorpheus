using Easy.Common.Extensions;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Readers;
using ScottPlot;
using SharpLearning.Containers.Extensions;
using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Drawing;
using System.IO;
using System.Linq;
using EngineLayer.CombinatorialSearch;
using TopDownProteomics;
using TorchSharp;
using UsefulProteomicsDatabases;
using IsotopicEnvelope = MassSpectrometry.IsotopicEnvelope;

namespace TaskLayer
{
    public class MMGPTMD
    {
        public Dictionary<int, Modification> ModsDictionary { set; get; }
        public List<Product> TheoreticalProducts { set; get; }

        public List<Tuple<List<MatchedFragmentIon>, int>> SearchResults { set; get; }
        public List<Tuple<List<MatchedFragmentIon>, int>> FinalResult { set; get; }
        public List<MatchedFragmentIon> NoModSearch { get; set; }
        public List<Protein> ProteinDb { get; set; }
        public PeptideWithSetModifications PeptideToSearch { get; set; }
        public IsotopicEnvelope[] DeconvolutedScan { get; set; }
        public List<MsDataScan> DataScans { get; set; }
        public List<MsDataScan> Ms1Scans { get; set; }
        public Modification[] UniModDb { get; set; }
        public List<PeptideWithSetModifications> PeptidesFromDb { get; set; }
        public Protein Peptide { get; set; }

        public MMGPTMD(int scanNumber, string filePath = @"D:\08-30-22_bottomup\test.mzML",
            string dbPath = @"D:\08-30-22_bottomup\database_example.fasta")
        {
            ModsDictionary = new Dictionary<int, Modification>();


            ProteinDb = ProteinDbLoader.LoadProteinFasta(dbPath,
                generateTargets: true,
                decoyType: DecoyType.None, isContaminant: false,
                out List<string> errors);

            var msDataFile = MsDataFileReader.GetDataFile(filePath);
            var msDataScans = msDataFile.LoadAllStaticData();

            var ms2Scans =
                from scans in msDataScans.Scans
                where scans.MsnOrder == 2
                select scans;

            DataScans = ms2Scans.ToList();

            var ms1Scans =
                from scans in msDataScans.Scans
                where scans.MsnOrder == 1
                select scans;

            Ms1Scans = ms1Scans.ToList();

            var deconvoluted = new Deconvoluter(DeconvolutionType.ClassicDeconvolution,
                new ClassicDeconvolutionParameters(1, 6, 10, 5));

            DeconvolutedScan = deconvoluted.Deconvolute(Ms1Scans[scanNumber]).ToArray();

            var tolerance = new PpmTolerance(100);
            PeptidesFromDb = new();
            for (int i = 0; i < ProteinDb.Count; i++)
            {
                var candidate = new PeptideWithSetModifications(
                    protein: ProteinDb[i],
                    new DigestionParams(), oneBasedStartResidueInProtein: 1,
                    oneBasedEndResidueInProtein: ProteinDb[i].BaseSequence.Length,
                    cleavageSpecificity: CleavageSpecificity.Full,
                    peptideDescription: String.Empty, missedCleavages: 1,
                    allModsOneIsNterminus: new Dictionary<int, Modification>(),
                    numFixedMods: 0);

                PeptidesFromDb.Add(candidate);
            }

            SearchResults = new List<Tuple<List<MatchedFragmentIon>, int>>();
            foreach (var pepide in ProteinDb)
            {
                ModsDictionary = new Dictionary<int, Modification>();

                Peptide = pepide;

                CommonFixedPTMs();

                GetPeptideWithSetModifications();

                TheoreticalProducts = new List<Product>();

                PeptideToSearch.Fragment(dissociationType: DissociationType.HCD,
                    fragmentationTerminus: FragmentationTerminus.Both,
                    products: TheoreticalProducts);

                NoModSearch = MetaMorpheusEngine.MatchFragmentIons(new Ms2ScanWithSpecificMass(DataScans[scanNumber],
                        DataScans[scanNumber].SelectedIonMZ.Value,
                        DataScans[scanNumber].SelectedIonChargeStateGuess.Value, filePath, new CommonParameters()),
                    TheoreticalProducts, new CommonParameters());


                SearchResults.Add(
                    new Tuple<List<MatchedFragmentIon>, int>(item1: NoModSearch, item2: NoModSearch.Count()));

                foreach (var search in SearchResults)
                {
                    Console.WriteLine("Experimental: " + search.Item1.Count + "| Theoretical: " +
                                      TheoreticalProducts.Count + "| Candidate: ");
                }

                Console.WriteLine("----------------");
            }

            UniModDb = UsefulProteomicsDatabases.Loaders.LoadUnimod(@"unimod.xml").ToArray();


        }

        public void GetPeptideWithSetModifications()
        {
            PeptideToSearch = new PeptideWithSetModifications(
                protein: Peptide,
                new DigestionParams(), oneBasedStartResidueInProtein: 1,
                oneBasedEndResidueInProtein: Peptide.BaseSequence.Length,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: String.Empty, missedCleavages: 1,
                allModsOneIsNterminus: ModsDictionary,
                numFixedMods: 0);
        }

        public void CommonFixedPTMs()
        {
            var peptideSequence = Peptide.BaseSequence.ToCharArray();

            int counter = 1;
            int PhosphoTwo = 0;
            for (int i = 0; i < peptideSequence.Length; i++)
            {
                if (peptideSequence[i].Equals('C'))
                {
                    ModsDictionary.Add(counter,
                        new Modification(_originalId: "Carbamidomethyl", _monoisotopicMass: 57.021464));
                    counter++;
                    ModsDictionary.Add(counter,
                        new Modification(_originalId: "Deamidation", _monoisotopicMass: -17.026549));
                    //ModsDictionary.Add(counter,
                    //    new Modification(_originalId: "Oxidation", _monoisotopicMass: 15.994915));
                    counter++;
                }
                //else if (peptideSequence[i].Equals('M'))
                //{
                //    ModsDictionary.Add(counter,
                //        new Modification(_originalId: "Acetylation", _monoisotopicMass: 42.010565));
                //    counter++;
                //    ModsDictionary.Add(counter,
                //        new Modification(_originalId: "Oxidation", _monoisotopicMass: 15.994915));
                //    counter++;
                //}
                else if (peptideSequence[i].Equals('S') && PhosphoTwo < 2)
                {
                    PhosphoTwo++;
                    counter++;
                }
                else if (peptideSequence[i].Equals('S') && PhosphoTwo >= 2)
                {
                    ModsDictionary.Add(counter, new Modification(_originalId: "Phospho", _monoisotopicMass: 79.9799));
                    counter++;
                }
                else
                {
                    counter++;
                }
            }
        }



        public static void MatchSpectra(string filePath = @"D:\08-30-22_bottomup\test.mzML")
        {
            var msDataFile = MsDataFileReader.GetDataFile(filePath);
            var msDataScans = msDataFile.LoadAllStaticData();

            var ms2Scans =
                from scans in msDataScans.Scans
                where scans.MsnOrder == 2
                select scans;

            var ms2ScansList = ms2Scans.ToList();
            var ms1Scans =
                from scans in msDataScans.Scans
                where scans.MsnOrder == 1
                select scans;

            var phosphoMod = new Mods().AAsMonoIsotopic;

            Dictionary<int, Modification> modsDictionary = new();
            Dictionary<int, Modification> acetylation = new();
            Dictionary<int, Modification> oxidation = new();
            Dictionary<int, Modification> oxiAndAcetyl = new();
            oxiAndAcetyl.Add(1, new Modification(_monoisotopicMass: 42.010565));
            oxiAndAcetyl.Add(2, new Modification(_monoisotopicMass: 15.994915));
            //modsDictionary.Add(1, new Modification(_monoisotopicMass: 57.021464));
            //modsDictionary.Add(23, new Modification(_monoisotopicMass:79.99));
            acetylation.Add(1, new Modification(_monoisotopicMass: 42.010565)); //acetylation
            oxidation.Add(1, new Modification(_monoisotopicMass: 15.994915)); //oxidation
            //GlobalVariables.AllModsKnown.Select(x => new KeyValuePair<int, Modification>(1, x));

            //get proteins from fasta database
            var proteinDB = ProteinDbLoader.LoadProteinFasta(@"D:\08-30-22_bottomup\database_example.fasta",
                generateTargets: true,
                decoyType: DecoyType.None, isContaminant: false, out List<string> errors);

            var mod = new PeptideWithSetModifications(
                protein: proteinDB[21],
                new DigestionParams(), oneBasedStartResidueInProtein: 1,
                oneBasedEndResidueInProtein: proteinDB[21].BaseSequence.Length,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: String.Empty, missedCleavages: 1,
                allModsOneIsNterminus: oxiAndAcetyl,
                numFixedMods: 0);

            var products = new List<Product>();

            mod.Fragment(dissociationType: DissociationType.HCD, fragmentationTerminus: FragmentationTerminus.Both,
                products: products);

            var matched = MetaMorpheusEngine.MatchFragmentIons(new Ms2ScanWithSpecificMass(ms2ScansList[21],
                    ms2ScansList[21].SelectedIonMZ.Value,
                    ms2ScansList[21].SelectedIonChargeStateGuess.Value, filePath, new CommonParameters()),
                products, new CommonParameters());

            List<Product> reducedProducts = new List<Product>(products);

            foreach (var product in products)
            {
                foreach (var match in matched)
                {
                    if (product.Annotation.Equals(match.NeutralTheoreticalProduct.Annotation))
                    {
                        reducedProducts.Remove(product);
                    }
                }
            }
        }

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

        public static Dictionary<string, double> GetModsDictionary()
        {
            var mods = GetModsFromGptmdThing();
            var aaDict = new Mods();

            var modsDictonary = new Dictionary<string, double>();

            foreach (var mod in mods)
            {
                var aminoacid = aaDict.AAsMonoIsotopic.TryGetValue(mod.Target.ToString(), out double residueMass);
                if (residueMass > 0 && modsDictonary.ContainsValue((double)mod.MonoisotopicMass + residueMass) == false)
                {
                    modsDictonary.Add(key: mod.Target.ToString() + " " + mod.OriginalId,
                        value: (double)mod.MonoisotopicMass + residueMass);
                }
            }

            return modsDictonary;
        }

        public static Dictionary<string, double> GetModsDictionaryNoAA()
        {
            var mods = GetModsFromGptmdThing();
            var aaDict = new Mods();

            var modsDictonary = new Dictionary<string, double>();

            foreach (var mod in mods)
            {

                modsDictonary.Add(key: mod.OriginalId + " on " + mod.Target, value: (double)mod.MonoisotopicMass);

            }

            return modsDictonary;
        }

        public static Dictionary<string, double> GetCandidates()
        {
            var mods = GetModsDictionary();
            var aa = new Mods().AAsMonoIsotopic;

            List<Dictionary<string, double>> candidateList = new();

            candidateList.Add(aa);
            candidateList.Add(mods);

            Dictionary<string, double> candidates = candidateList.SelectMany(dict => dict)
                .ToDictionary(pair => pair.Key, pair => pair.Value);

            return candidates;
        }

        public static List<Mods> MultiModDiscovery(string filePath)
        {

            var msDataFile = MsDataFileReader.GetDataFile(filePath);
            var msDataScans = msDataFile.LoadAllStaticData();

            var ms2Scans =
                from scans in msDataScans.Scans
                where scans.MsnOrder == 2
                select scans;

            var ms1Scans =
                from scans in msDataScans.Scans
                where scans.MsnOrder == 1
                select scans;

            //var mods = new Mods(@"Data\unimod.xml");


            var mods2 = GlobalVariables.AllModsKnownDictionary;
            foreach (var scan in ms2Scans)
            {
                double[] mz = scan.MassSpectrum.XArray;
                //double[] deltaMz = new double[] { mz.Length };

                //torch.Tensor mz =  torch.tensor(scan.MassSpectrum.XArray);
                List<List<double>> deltaMz = new();
                //List < torch.Tensor > deltaMzTensors = new();

                for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                {
                    List<double> tempDeltaMzList = new();

                    for (int k = i + 1; k < scan.MassSpectrum.XArray.Length - 1; k++)
                    {
                        var tempDeltaMz = Math.Abs(mz[i] - mz[k]);
                        tempDeltaMzList.Add(tempDeltaMz);
                    }

                    deltaMz.Add(tempDeltaMzList);
                }

                var mods = GetCandidates();
                var tolerance = new PpmTolerance(20); //Tolerance is better than more or less 
                List<List<string>> seq = new();
                List<List<KeyValuePair<string, double>>> possibleSequences = new();
                foreach (var delta in deltaMz)
                {
                    List<string> candidateSeq = new List<string>();

                    var possibleSeq =
                        from i in delta
                        where i > 40
                        select i;
                    List<KeyValuePair<string, double>> possibleMod = new();

                    foreach (var mod in possibleSeq)
                    {

                        foreach (var candidate in mods)
                        {
                            if (tolerance.Within(mod, candidate.Value))
                            {
                                possibleMod.Add(candidate);
                            }
                        }
                    }

                    possibleSequences.Add(possibleMod);
                    seq.Add(candidateSeq);
                }


                foreach (var sequence in possibleSequences)
                {
                    if (sequence.Count > 1)
                    {
                        Console.Write("[");
                        foreach (var candidate in sequence)
                        {
                            Console.Write(candidate.Key + " , ");
                        }

                        Console.Write("]");
                    }
                    else
                    {
                        foreach (var candidate in sequence)
                        {
                            Console.Write(candidate.Key + " | ");
                        }
                    }

                    Console.WriteLine();
                }

                ////foreach (var delta in deltaMz)
                ////{
                ////    foreach (var residue in delta)
                ////    {
                ////        var AADict = new Mods().AAsMonoIsotopic.ToList();
                ////        var modList = mods.ToImmutableList();
                ////        var roundedResidue = (double)Math.Round(residue, 2);



                ////        for (int i = 0; i < AADict.Count(); i++)
                ////        {
                ////            if (tolerance.Within(residue, AADict[i].Value))
                ////            {
                ////                Console.Write(AADict[i].Key + " | ");
                ////                break;
                ////            }
                ////        }
                ////        for (int j = 0; j < modList.Count(); j++)
                ////        {
                ////            if (tolerance.Within(residue, modList[j].Key))
                ////            {
                ////                Console.Write(modList[j].Value + " | ");
                ////                break;
                ////            }
                ////        }   


                ////        //if(tolerance.Within(residue, AADict.Get))
                ////        ////Residue.GetResidue().MonoisotopicMass; GETS THE MONOISOTOPIC
                ////        //foreach (var AA in AADict)
                ////        //{
                ////        //    if(tolerance.Within(residue, AA.Value))
                ////        //    {
                ////        //        Console.Write(AA.Key + " | ");
                ////        //    }
                ////        //    else
                ////        //    {
                ////        //        var aaResName = mods.Values.ToArray();
                ////        //        var aaMassArray= mods.Keys.Select(x => Math.Round(x, 3)).ToArray();
                ////        //        var binarySearch = Array.BinarySearch(aaMassArray, Math.Round(residue, 3));
                ////        //        if (binarySearch >= 0)
                ////        //        {
                ////        //            Console.WriteLine(aaResName[binarySearch]);
                ////        //        }
                ////        //    }
                ////        //}
                ////    }
                break;
            }

            return new List<Mods>();
        }

        public static void UpdateTheFilteredPsmFile(MsDataFile msDataFile, string psmFilePath)
        {
            List<FilteredPsmTSV> psms = ReadFilteredPsmTSVShort(psmFilePath);

            MsDataScan[] dataScans = msDataFile.GetMsDataScans();

            int precursorNumber = 1;
            int scanNumber = 2;
            //Use list instead of Scans GEScansList()
            for (int i = 0; i < psms.Count(); i++)
            {
                psms[i].PrecursorScanNumber = (precursorNumber).ToString();
                psms[i].ScanNumber = (scanNumber).ToString();
                precursorNumber = precursorNumber + 2;
                scanNumber = scanNumber + 2;
            }

            WriteUpdatedFilteredPsmsToTSV(psms, psmFilePath);
        }

        public static void WriteUpdatedFilteredPsmsToTSV(List<FilteredPsmTSV> filteredPsms, string path)
        {
            using (var writer = new StreamWriter(path))
            {
                //makes this an enum?
                string header = "File Name\tScan Number\tPrecursor Scan Number\tScore\tBase Sequence\t" +
                                "Full Sequence\tMods\tMods Count\tProtein Accession\t " +
                                "Protein Name\tGene Name\tOrganism Name\t" +
                                "Start and End Residues in Protein\t" +
                                "Matched Ion Series\tMatched Ion Counts\tPrecursor Mass\tCharge";

                writer.WriteLine(header);
                foreach (var psm in filteredPsms)
                {
                    string[] row = new[]
                    {
                        psm.FileName,
                        psm.ScanNumber,
                        psm.PrecursorScanNumber,
                        psm.Score,
                        psm.BaseSeq,
                        psm.FullSeq,
                        psm.Mods,
                        psm.ModsCount,
                        psm.ProteinAccession,
                        psm.ProteinName,
                        psm.GeneName,
                        psm.OrganismName,
                        psm.StartAndEndResiduesInProtein,
                        psm.MatchedIonSeries,
                        psm.MatchedIonCounts,
                        psm.PrecursorMass,
                        psm.Charge
                    };
                    writer.WriteLine(string.Join('\t', row));
                }
            }
        }

        public static void CreateMzML(MsDataScan[] scans, SourceFile sourceFile, string path)
        {
            SourceFile sourceFile1 = new("no nativeID format", "mzML format",
                null, null, filePath: @"K:\08-30-22_bottomup\fractionated_search\Task3-SearchTask\file_example.mzML",
                null);

            MsDataFile genericFile = new GenericMsDataFile(scans: scans, sourceFile: sourceFile);

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(genericFile, path, true);
        }

        public static (MsDataScan[], MsDataFile) ExtractScansAndSourceFile(List<FilteredPsmTSV> psms,
            List<string> filePaths)
        {
            List<MsDataFile> loadedFiles = new();

            foreach (var file in filePaths)
            {
                loadedFiles.Add(Readers.MsDataFileReader.GetDataFile(file).LoadAllStaticData());
            }

            List<string> fileName = new();

            foreach (var name in loadedFiles)
            {
                fileName.Add(name.FilePath);
            }

            List<Tuple<string, MsDataFile>> tupleList = new();

            MsDataFile dataFile = loadedFiles[0];

            for (int i = 0; i < loadedFiles.Count(); i++)
            {
                tupleList.Add(new Tuple<string, MsDataFile>(item1: fileName[i], item2: loadedFiles[i]));
            }

            var dict = tupleList.ToImmutableDictionary(x => x.Item1, x => x.Item2);

            List<MsDataScan> scanList = new List<MsDataScan>();
            int counter = 1;
            int precursorNumber = 1;
            double retentionTime = 1;
            double injectionTime = 1;
            foreach (var psm in psms)
            {
                foreach (var file in dict)
                {
                    if (file.Key.Contains(psm.FileName))
                    {
                        //var ms1 = file.Value.GetMS1Scans().First();
                        MsDataScan ms2 = file.Value.GetOneBasedScan(int.Parse(psm.ScanNumber));
                        var ms1 = file.Value.GetOneBasedScan(int.Parse(psm.PrecursorScanNumber));
                        double[,] mzIntensitiesMS1 = new double[2, ms1.MassSpectrum.XArray.Length];
                        double[,] mzIntensitiesMS2 = new double[2, ms2.MassSpectrum.XArray.Length];

                        // Two 1-D array to One 2-D array
                        for (int i = 0; i < ms1.MassSpectrum.XArray.Length; i++)
                        {
                            mzIntensitiesMS1[0, i] = ms1.MassSpectrum.XArray[i];
                            mzIntensitiesMS1[1, i] = ms1.MassSpectrum.YArray[i];
                        }

                        for (int i = 0; i < ms2.MassSpectrum.XArray.Length; i++)
                        {
                            mzIntensitiesMS2[0, i] = ms2.MassSpectrum.XArray[i];
                            mzIntensitiesMS2[1, i] = ms2.MassSpectrum.YArray[i];
                        }

                        MzSpectrum spectrumMS1 = new MzSpectrum(mzIntensitiesMS1);
                        MzSpectrum spectrumMS2 = new MzSpectrum(mzIntensitiesMS2);
                        //Recreate Scan
                        MsDataScan scanMS1 = new MsDataScan(massSpectrum: spectrumMS1, oneBasedScanNumber: counter,
                            msnOrder: ms1.MsnOrder, isCentroid: ms1.IsCentroid, polarity: ms1.Polarity,
                            retentionTime: retentionTime, scanWindowRange: ms1.ScanWindowRange,
                            scanFilter: ms1.ScanFilter, mzAnalyzer: ms1.MzAnalyzer,
                            totalIonCurrent: ms1.TotalIonCurrent, injectionTime: injectionTime,
                            noiseData: ms1.NoiseData,
                            nativeId: "controllerType=0 controllerNumber=1 scan=" + counter.ToString(),
                            selectedIonMz: ms1.SelectedIonMZ,
                            selectedIonChargeStateGuess: ms1.SelectedIonChargeStateGuess,
                            selectedIonIntensity: ms1.SelectedIonIntensity,
                            isolationMZ: ms1.IsolationMz, isolationWidth: ms1.IsolationWidth,
                            dissociationType: ms1.DissociationType,
                            hcdEnergy: ms1.HcdEnergy);

                        counter = counter + 1;
                        retentionTime++;
                        injectionTime++;

                        MsDataScan scanMS2 = new MsDataScan(massSpectrum: spectrumMS2, oneBasedScanNumber: counter,
                            msnOrder: ms2.MsnOrder, isCentroid: ms2.IsCentroid, polarity: ms2.Polarity,
                            retentionTime: retentionTime, scanWindowRange: ms2.ScanWindowRange,
                            scanFilter: ms2.ScanFilter, mzAnalyzer: ms2.MzAnalyzer,
                            totalIonCurrent: ms2.TotalIonCurrent, injectionTime: injectionTime,
                            noiseData: ms2.NoiseData,
                            nativeId: "controllerType=0 controllerNumber=1 scan=" + counter.ToString(),
                            selectedIonMz: ms2.SelectedIonMZ,
                            selectedIonChargeStateGuess: ms2.SelectedIonChargeStateGuess,
                            selectedIonMonoisotopicGuessMz: ms2.SelectedIonMonoisotopicGuessMz,
                            selectedIonIntensity: ms2.SelectedIonIntensity,
                            isolationMZ: ms2.IsolationMz, isolationWidth: ms2.IsolationWidth,
                            dissociationType: ms2.DissociationType, oneBasedPrecursorScanNumber: counter - 1,
                            hcdEnergy: ms2.HcdEnergy);


                        counter = counter + 1;
                        precursorNumber++;
                        retentionTime++;
                        injectionTime++;
                        //tempFile.SetOneBasedScanNumber(counter);

                        scanList.Add(scanMS1);
                        scanList.Add(scanMS2);
                        break;
                    }
                }
            }

            MsDataScan[] scansArray = Enumerable.Select(scanList, x => x).ToArray();

            return (scansArray, dataFile);
        }

        public static void WriteFastaDBFromFilteredPsm(List<FilteredPsmTSV> psms, string path)
        {
            List<Protein> proteinList = new List<Protein>();

            foreach (var protein in psms)
            {
                proteinList.Add(new Protein(
                    sequence: protein.BaseSeq,
                    accession: protein.ProteinAccession,
                    name: protein.ProteinName,
                    organism: protein.OrganismName
                ));
            }

            ProteinDbWriter.WriteFastaDatabase(proteinList, path, "");
        }


        public static IEnumerable<PsmFromTsv> FilterPsm(List<PsmFromTsv> psms)
        {
            IEnumerable<PsmFromTsv> filteredPsms =
                from psm in psms
                where psm.QValue <= 0.0001 && psm.PEP <= 0.0001 && psm.DecoyContamTarget.Equals("T") //&&
                      //PsmFromTsv.ParseModifications(psm.FullSequence).Count() >= 2
                //where psm.QValue <= 0.001 && psm.PEP <= 0.001 && psm.DecoyContamTarget.Equals("T") &&
                //      PsmFromTsv.ParseModifications(psm.FullSequence).Count() >= 2
                select psm;

            return filteredPsms.Take(10).Distinct();
        }

        public static void WriteFilteredPsmsToTSVFull(List<PsmFromTsv> filteredPsms, string path)
        {
            using (var writer = new StreamWriter(path))
            {
                //makes this an enum?
                string header =
                    "File Name\tScan Number\tScan Retention Time\tNum Experimental Peaks\tTotal Ion Current\t" +
                    "Precursor Scan Number\tPrecursor Charge\tPrecursor MZ\tPrecursor Mass\tScore\tDelta Score\t" +
                    "Notch\tBase Sequence\tFull Sequence\tEssential Sequence\tAmbiguity Level\tPSM Count (unambiguous, <0.01 q-value)\t" +
                    "Mods\tMods Chemical Formulas\tNum Variable Mods\tMissed Cleavages\tPeptide Monoisotopic Mass\tMass Diff (Da)\t" +
                    "Mass Diff (ppm)\tProtein Accession\tProtein Name\tGene Name\tOrganism Name\tIntersecting Sequence Variations\t" +
                    "Identified Sequence Variations\tSplice Sites\tContaminants\tDecoy\tPeptide Description\tStart and End Residues In Protein\t" +
                    "Previous Amino Acid\tNext Amino Acid\tTheoreticals Searched\tDecoy/Contaminant/Target\tMatched Ion Series\tMatched Ion Mass-To-Charge Ratios\t" +
                    "Matched Ion Mass Diff (Da)\tMatched Ion Mass Diff (Ppm)\tMatched Ion Intensities\tMatched Ion Counts\tLocalized Scores\tImprovement Possible\t" +
                    "Cumulative Target\tCumulative Decoy\tCumulative Target Notch\tCumulative Decoy Notch\tQValue\tQValue Notch\tPEP\tPEP_QValue";

                writer.WriteLine(header);
                foreach (var psm in filteredPsms)
                {
                    string[] row = new[]
                    {
                        String.Join('-', psm.FileNameWithoutExtension.Split('-').SkipLast(1)),
                        psm.Ms2ScanNumber.ToString(),
                        psm.PrecursorScanNum.ToString(),
                        psm.RetentionTime.ToString(),
                        string.Join(' ', Enumerable.Select(psm.MatchedIons, x => x.Annotation).ToArray()),
                        psm.MatchedIons.Count().ToString()
                    };
                    writer.WriteLine(string.Join('\t', row));
                }
            }
        }

        public static void WriteFilteredPsmsToTSV(List<PsmFromTsv> filteredPsms, string path)
        {
            using (var writer = new StreamWriter(path))
            {
                //makes this an enum?
                string header = "File Name\tScan Number\tPrecursor Scan Number\tScore\tBase Sequence\t" +
                                "Full Sequence\tMods\tMods Count\tProtein Accession\t " +
                                "Protein Name\tGene Name\tOrganism Name\t" +
                                "Start and End Residues in Protein\t" +
                                "Matched Ion Series\tMatched Ion Counts\tPrecursor Mass\tCharge";

                writer.WriteLine(header);
                foreach (var psm in filteredPsms)
                {
                    string[] row = new[]
                    {
                        String.Join('_', psm.FileNameWithoutExtension.Split(new [] {'_'})),
                        psm.Ms2ScanNumber.ToString(),
                        psm.PrecursorScanNum.ToString(),
                        psm.Score.ToString(),
                        psm.BaseSeq,
                        psm.FullSequence,
                        String.Join(", ",
                            PsmFromTsv.ParseModifications(psm.FullSequence).Select(x =>
                                    string.Concat(string.Join(", ", x.Value.ToArray()), ", ",
                                        Math.Abs(x.Key).ToString()))
                                .ToArray()),
                        PsmFromTsv.ParseModifications(psm.FullSequence).Count().ToString(),
                        psm.ProteinAccession,
                        psm.ProteinName,
                        psm.GeneName,
                        psm.OrganismName,
                        psm.StartAndEndResiduesInProtein,
                        string.Join(' ', Enumerable.Select(psm.MatchedIons, x => x.Annotation).ToArray()),
                        psm.MatchedIons.Count().ToString(),
                        psm.PrecursorMass.ToString(),
                        psm.PrecursorCharge.ToString(),
                    };
                    writer.WriteLine(string.Join('\t', row));
                }
            }
        }

        public static List<PsmFromTsv> ReadPSMTSVFull(string psmPath)
        {
            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmPath, out var warnings);
            return parsedPsms;

        }

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

    public class ScalableModSearch
    {
        public List<Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>> CombinationOfModificationsMap { get; set; }
        public List<torch.Tensor> ComboMassTensors { get; set; }

        public ScalableModSearch(MsDataFile msDataFile, List<FilteredPsmTSV> psmList)
        {
            SetCandidateList(msDataFile, psmList);
        }

        public List<torch.Tensor> SetComboMassTensors()
        {
            foreach (var peptide in CombinationOfModificationsMap)
            {
                var peptideComboMassList = peptide.Select(x => Enumerable.Select(x.Value, ion => ion.Mz).ToArray()).ToList();
                //var peptideComboMassTensor = torch.tensor(peptideComboMassList.Select((x => x)).ToArray());
                foreach (var mod in peptide)
                {

                }
            }

            return new List<torch.Tensor>();
        }

        public void SetCandidateList(MsDataFile msDataFile, List<FilteredPsmTSV> psmList)
        {
            CombinationOfModificationsMap = GetAllComboMods(msDataFile, psmList);
        }


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

        public static void GetModsPresentInPsm(string psmPath, out List<Dictionary<int, Modification>> modsFromPsm)
        {
            var psms = ReadFilteredPsmTSVShort(psmPath);
            modsFromPsm = new List<Dictionary<int, Modification>>();
            var unimod = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation: "unimod.xml").ToList();
            foreach (var psm in psms)
            {

                var dictionary = MMGPTMD.GetModsDictionaryNoAA();
                var modsDict = new Dictionary<int, Modification>();
                var psmSplit = psm.Mods.Split(", ");
                for (int i = 0; i < psmSplit.Length; i++)
                {
                    if (psmSplit[i].Contains("Common Artifact:"))
                        psmSplit[i] = psmSplit[i].Replace("Common Artifact:", "");
                    else if (psmSplit[i].Contains("Common Fixed:"))
                        psmSplit[i] = psmSplit[i].Replace("Common Fixed:", "");
                    else if (psmSplit[i].Contains("Common Biological:"))
                        psmSplit[i] = psmSplit[i].Replace("Common Biological:", "");
                }

                for (int i = 0; i < psmSplit.Length / 2; i = i + 2)
                {

                    double monoIsotopicMass = 0;
                    char target = psmSplit[i][psmSplit[i].Length - 1];
                    dictionary.TryGetValue(psmSplit[i], out monoIsotopicMass);
                    modsDict.Add(Math.Abs(int.Parse(psmSplit[i + 1]) + 1),
                        new Modification(_monoisotopicMass: monoIsotopicMass));
                }

                modsFromPsm.Add(modsDict);
            }
        }

        public void GetTopScoreAndSavePng(
            List<Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>> getAllComboModsResults,
            string pathToSavePlots)
        {
            var msDataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML");

            List<KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>> topCount = new();

            foreach (var result in getAllComboModsResults)
            {
                foreach (var pair in result.OrderByDescending(x => x.Value.Count))
                {
                    topCount.Add(pair);
                    break;
                }
            }
            //ScalableModSearch.GetModsPresentInPsm(@"D:\08-30-22_bottomup\example.psmtsv", out testing);

            foreach (var result in topCount)
            {
                MsDataScan scan = msDataFile.GetOneBasedScan(int.Parse(result.Key.PeptideDescription));
                double[] mz = scan.MassSpectrum.XArray;
                double[] intesities = Enumerable.Select(scan.MassSpectrum.YArray, x => ScalableModSearch.NormalizeArray(
                    scan.MassSpectrum.YArray, x)).ToArray();

                var plt = new ScottPlot.Plot(2000, 800);
                for (int i = 0; i < mz.Length; i++)
                {
                    var vlines = new ScottPlot.Plottable.VLineVector();
                    vlines.Xs = new[] { mz[i] };
                    vlines.Max = intesities[i];
                    vlines.Color = Color.DarkGray;
                    //vlines.PositionLabel = true;
                    //vlines.PositionLabelBackground = vlines.Color;
                    plt.Add(vlines);

                }

                var rand = new Random();

                for (int i = 0; i < result.Value.Count(); i++)
                {
                    //var arrow = plt.AddArrow(result.Value[i].Mz, intesities[i] + 10, result.Value[i].Mz,
                    //    intesities[i] + randomValue, lineWidth: 0.5f);
                    //arrow.Label = result.Value[i].Annotation;
                    //arrow.Color = Color.Black;
                    if (result.Value[i].Annotation.ToCharArray()[0].Equals('y'))
                    {
                        var randomValue = rand.NextInt64((long)intesities[i] + 10, (long)intesities[i] + 25);
                        var label = plt.AddText(result.Value[i].Annotation,
                            result.Value[i].Mz, randomValue);// + randomValue+2);
                        var line = plt.AddVerticalLine(result.Value[i].Mz, Color.LightCoral);
                        line.Max = randomValue;
                        label.Color = Color.Red;
                    }
                    else
                    {
                        var randomValue = rand.NextInt64((long)intesities[i] + 20, (long)intesities[i] + 30);
                        var label = plt.AddText(result.Value[i].Annotation,
                            result.Value[i].Mz, randomValue);// + randomValue+2);
                        var line = plt.AddVerticalLine(result.Value[i].Mz, Color.LightBlue);
                        line.Max = randomValue;
                        label.Color = Color.Blue;
                    }
                }
                //foreach (var ion in result.Value)
                //{
                //    var arrow = plt.AddArrow(ion.Mz, ion.Intensity, ion.Mz, ion.Intensity);
                //    arrow.Label = ion.Annotation;
                //    arrow.Color = Color.Black;

                //    var label = plt.AddText(ion.Annotation, ion.Mz, ion.Intensity);
                //    label.Color = Color.Black;
                //    //plt.AddText(ion.Annotation, ion.Mz, ion.Intensity, color: Color.Red, size: 16);
                //    //plt.AddMarker(ion.Mz, ion.Intensity, MarkerShape.verticalBar,
                //    //    label: ion.Annotation);
                //}
                plt.AddAnnotation("Theoretical: " + result.Key.MonoisotopicMass + " Da\n" + $"Precursor * charge: {result.Key.Protein.Name} Da",
                    alignment: Alignment.UpperRight);
                plt.AddAnnotation(result.Key.Protein.Accession, Alignment.UpperLeft);
                plt.Title(result.Key.FullSequence + ", " + result.Value.Count() + " matches", size: 10);
                plt.YLabel("Intensity");
                plt.XLabel("m/Z");
                plt.SetAxisLimitsY(0, intesities.Max() + 5);
                plt.XAxis.Layout(padding: 20);
                plt.YAxis.Layout(padding: 20);
                //plt.YAxis.LabelStyle(fontSize: 24);
                //plt.XAxis.LabelStyle(fontSize: 24);
                //plt.XAxis2.LabelStyle(fontSize: 36);

                //plt.YAxis.TickLabelStyle(fontSize: 18);
                //plt.XAxis.TickLabelStyle(fontSize: 18);

                //plt.YAxis.MajorGrid(lineWidth: 2);
                //plt.XAxis.MajorGrid(lineWidth: 2);
                plt.SaveFig(@$"{pathToSavePlots}{result.Key.BaseSequence}.png", scale: 10);

            }
        }

        /// <summary>
        /// Nic's talk on aug1
        /// </summary>
        /// <param name="dataFile"></param>
        /// <param name="filteredPsmList"></param>
        /// <returns></returns>
        public static List<Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>> NewGetAllComboMods(
            MsDataFile dataFile,
            List<FilteredPsmTSV> filteredPsmList)
        {
            List<Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>> finalCandidates = new();
            var matches = new List<Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>>();

            foreach (var psm in filteredPsmList)
            {
                char[] baseSequence = psm.BaseSeq.ToCharArray();
                var spectrum = dataFile.GetOneBasedScan(int.Parse(psm.ScanNumber));
                var precursorSelectedPrecursorMass =
                    ((spectrum.SelectedIonMonoisotopicGuessMz) * spectrum.SelectedIonChargeStateGuess); // -
                //(spectrum.SelectedIonChargeStateGuess * 1.00727647);

                //todo make peptide from psm Base seq and digest with top down protease, that way I have all the properties of the object. Treating the peptide as a protein

                var mods =
                    Loaders.LoadUnimod(
                            @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                        .ToList();
                var commonModsFromToml = MMGPTMD.GetModsFromGptmdThing().ToList();
                commonModsFromToml.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")));

                var peptideProteinDigest =
                    new Protein(psm.BaseSeq, psm.ProteinAccession).Digest(new DigestionParams(protease: "top-down"),
                        allKnownFixedModifications: new List<Modification>()
                            { mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")) },
                        variableModifications: new List<Modification>());



                var modifiedPeptide = new PeptideWithSetModifications(
                   peptideProteinDigest.First().Protein, new DigestionParams(), oneBasedStartResidueInProtein: 1,
                    oneBasedEndResidueInProtein: peptideProteinDigest.First().BaseSequence.Length,
                    CleavageSpecificity.Full, peptideDescription: spectrum.OneBasedScanNumber.ToString(), missedCleavages: 1,
                    allModsOneIsNterminus: peptideProteinDigest.First().AllModsOneIsNterminus, numFixedMods: peptideProteinDigest.First().NumFixedMods);

                var nakedProteinMass = new Protein(psm.BaseSeq, psm.ProteinAccession).Digest(new DigestionParams(protease: "top-down"),
                    allKnownFixedModifications: new List<Modification>(),
                    variableModifications: new List<Modification>());

                double TheoPrecursorMassWithMods = modifiedPeptide.MonoisotopicMass;

                var products = new List<Product>();

                modifiedPeptide.Fragment(dissociationType: DissociationType.HCD,
                    fragmentationTerminus: FragmentationTerminus.Both, products: products);

                var match = MetaMorpheusEngine.MatchFragmentIons(new Ms2ScanWithSpecificMass(
                    mzLibScan: spectrum, precursorMonoisotopicPeakMz: spectrum.SelectedIonMZ.Value,
                    precursorCharge: spectrum.SelectedIonChargeStateGuess.Value, fullFilePath: @"D:\08-30-22_bottomup\test.mzML\",
                    new CommonParameters()), products, new CommonParameters());

                var deltaMass = precursorSelectedPrecursorMass - modifiedPeptide.MonoisotopicMass;

                int numberOfMaxMods = 3;

                int numberOfAttempts = 15;

                var tolerance = new PpmTolerance(15);

                List<Modification> fixedMod = new List<Modification>()
                    { mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")) };

                var possibles = from one in mods
                                from two in mods
                                //from three in commonModsFromToml

                                where tolerance.Within((double)deltaMass, (one.MonoisotopicMass.Value + two.MonoisotopicMass.Value))// + three.MonoisotopicMass.Value))
                                && nakedProteinMass.First().BaseSequence.Contains(one.Target.ToString()) && nakedProteinMass.First().BaseSequence.Contains(two.Target.ToString())
                                && !one.IdWithMotif.Equals(fixedMod.Select(x => x)) && !two.IdWithMotif.Equals(fixedMod.Select(x => x))
                                //&& nakedProteinMass.First().BaseSequence.Contains(three.Target.ToString())

                                select new KeyValuePair<double, Modification[]>((one.MonoisotopicMass.Value+two.MonoisotopicMass.Value),//+three.MonoisotopicMass.Value),
                                        value: new[]{one,two});

                var groupedByMass = from mod in commonModsFromToml
                                                group mod by mod.MonoisotopicMass
                                                into newGroup
                                                select newGroup;

                var groupedModsByOriginalId = from one in commonModsFromToml
                    group one by one.OriginalId
                    into newGroup
                    select newGroup;
                    //select new KeyValuePair<double, Modification[]>((one.MonoisotopicMass.Value + two.MonoisotopicMass.Value),//+three.MonoisotopicMass.Value),
                    //    value: new[] { one, two });

                    var combos = (from grouped1 in groupedModsByOriginalId
                        from grouped2 in groupedModsByOriginalId
                        select new KeyValuePair<double, Modification[]>(
                            (grouped1.First().MonoisotopicMass.Value + grouped2.First().MonoisotopicMass.Value),
                            new[] { grouped1.First(), grouped2.First() })).OrderBy(x => x.Key);

                    //var combosSorted = combos.ToImmutableSortedDictionary();

                    foreach (var possibility in possibles)
                {
                    var candidate = possibility.Value.Select(x => x).ToList();

                    var candidatePeptide  =
                        new Protein(psm.BaseSeq, psm.ProteinAccession).Digest(new DigestionParams(protease: "top-down"),
                            allKnownFixedModifications: fixedMod,
                            variableModifications: candidate);

                    var candidateFragmentation = new PeptideWithSetModifications(
                        candidatePeptide.First().Protein, new DigestionParams(), oneBasedStartResidueInProtein: 1,
                        oneBasedEndResidueInProtein: candidatePeptide.First().BaseSequence.Length,
                        CleavageSpecificity.Full, peptideDescription: spectrum.OneBasedScanNumber.ToString(), missedCleavages: 1,
                        allModsOneIsNterminus: candidatePeptide.First().AllModsOneIsNterminus, numFixedMods: candidatePeptide.First().NumFixedMods);

                    candidateFragmentation.Fragment(dissociationType: DissociationType.HCD,
                        fragmentationTerminus: FragmentationTerminus.Both, products: products);

                    var candidateMatch = MetaMorpheusEngine.MatchFragmentIons(new Ms2ScanWithSpecificMass(
                        mzLibScan: spectrum, precursorMonoisotopicPeakMz: spectrum.SelectedIonMZ.Value,
                        precursorCharge: spectrum.SelectedIonChargeStateGuess.Value, fullFilePath: @"D:\08-30-22_bottomup\test.mzML\",
                        new CommonParameters()), products, new CommonParameters());
                    
                    var candidateDict = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();

                    candidateDict.Add(candidateFragmentation, candidateMatch);
                    matches.Add(candidateDict);

                }

                break;
            }
            return matches;

        }

        private static IEnumerable<IEnumerable<T>> CartesianProduct<T>(IEnumerable<IEnumerable<T>> sequences)
        {
            IEnumerable<IEnumerable<T>> emptyProduct = new[] { Enumerable.Empty<T>() };
            return sequences.Aggregate(
                emptyProduct,
                (accumulator, sequence) =>
                    from accseq in accumulator
                    from item in sequence
                    select accseq.Concat(new[] { item })
            );
        }

        public static List<Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>> GetAllComboMods(
            MsDataFile dataFile,
            List<FilteredPsmTSV> filteredPsmList)
        {
            List<Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>> finalCandidates = new();
            var modDB = MMGPTMD.GetModsDictionaryNoAA();
            modDB.Add("Carbamidomethyl on C", 57.021464);
            foreach (var psm in filteredPsmList)
            {
                char[] baseSequence = psm.BaseSeq.ToCharArray();
                var spectrum = dataFile.GetOneBasedScan(int.Parse(psm.ScanNumber));
                var precursorMass = ((spectrum.SelectedIonMonoisotopicGuessMz) * spectrum.SelectedIonChargeStateGuess);// -
                                                                                                                       //(spectrum.SelectedIonChargeStateGuess * 1.00727647);

                //todo make peptide from psm Base seq and digest with top down protease, that way I have all the properties of the object. Treating the peptide as a protein
                //from here

                var mods =
                    Loaders.LoadUnimod(
                        @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml").ToList();
                ////var mods = MMGPTMD.GetModsFromGptmdThing().ToList();
                var peptideProtein =
                    new Protein(psm.BaseSeq, psm.ProteinAccession).Digest(new DigestionParams(protease: "top-down"),
                        allKnownFixedModifications: new List<Modification>() { mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")) },
                        variableModifications: new List<Modification>());
                //to here 

                var peptide = new PeptideWithSetModifications(
                    protein: new Protein(sequence: psm.BaseSeq, accession: psm.ProteinAccession),
                    new DigestionParams(), oneBasedStartResidueInProtein: 1,
                    oneBasedEndResidueInProtein: psm.BaseSeq.Length, cleavageSpecificity: CleavageSpecificity.Full,
                    peptideDescription: spectrum.OneBasedScanNumber.ToString(),
                    missedCleavages: 1, allModsOneIsNterminus: new Dictionary<int, Modification>(),
                    numFixedMods: 0);

                double modMassCounter = 0;
                List<Dictionary<int, Modification>> candidateMods = new();

                Dictionary<int, Modification> preFilledKeyPairs = new Dictionary<int, Modification>();
                var deltaMass = precursorMass - peptide.MonoisotopicMass;

                for (int i = 1; i < (baseSequence.Length * 2) + 1; i++)
                {
                    preFilledKeyPairs.Add(i, new Modification());
                }

                //Carboamidomethylate all C's
                for (int i = 0; i < baseSequence.Length; i++)
                {
                    if (baseSequence[i].Equals('C'))
                    {
                        preFilledKeyPairs[i + 2] = new Modification(_originalId: "Carbamidomethylation on C",
                            _monoisotopicMass: 57.021464);
                        deltaMass = deltaMass - 57.021363;
                    }

                    //if first C, add ammonia loss
                    if (baseSequence[i].Equals('C') && i == 0)
                    {
                        preFilledKeyPairs[i + 1] = new Modification(_originalId: "Ammonia Loss on C",
                            _monoisotopicMass: -17.026549);
                        deltaMass = deltaMass - (-17.026549);
                    }
                }

                //var modsAvailable = from mod in modDB
                //                    where mod.Value.AlmostEqual((double)deltaMass, 2) 
                //                    select mod;

                //todo: try with the entire db

                var unimod = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation: "unimod.xml").ToList();


                var modsAvailable = from mod in modDB
                                    where baseSequence.Contains(mod.Key.ToCharArray()[mod.Key.ToCharArray().Length - 1])
                                    select mod;

                var cycles = true;
                var cycleNum = 0;
                //start adding the combos to the candidate list
                //while (cycles.Equals(true))
                //{
                var randomNumber = new Random();

                double originalDeltaMass = (double)deltaMass;

                var shuffledMods = modsAvailable.ToList();
                shuffledMods.Shuffle(randomNumber);
                foreach (var mod in shuffledMods)
                {
                    for (int i = 0; i < baseSequence.Length; i++)
                    {
                        if (baseSequence[i].Equals(mod.Key.ToCharArray()[mod.Key.ToCharArray().Length - 1]))
                        {
                            var newCandidate = new Dictionary<int, Modification>(preFilledKeyPairs);

                            newCandidate[i + 1] = new Modification(_originalId: mod.Key, _monoisotopicMass: mod.Value);
                            originalDeltaMass = originalDeltaMass + mod.Value;
                            //if (originalDeltaMass > 1500)
                            //    break;
                            //if (originalDeltaMass < cycleNum)
                            //{
                            candidateMods.Add(newCandidate);
                            //    cycles = false;
                            //    break;
                            //}
                        }

                        //if (cycles == false)
                        //    break;

                    }
                }
                cycleNum = cycleNum + 1;
                //}

                foreach (var mod in candidateMods)
                {
                    int numberOfMods = (mod.Keys.Count + 1);
                    for (int i = 1; i < numberOfMods; i++)
                    {
                        if (mod[i].OriginalId.IsNullOrEmpty())
                        {
                            mod.Remove(i);
                        }
                    }
                }

                Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> matches = new();

                //fragment and compare the ion count
                foreach (var candidate in candidateMods)
                {
                    //precursor mass embedded in the Protein's object name
                    var modifiedPeptide = new PeptideWithSetModifications(
                        protein: new Protein(sequence: peptide.BaseSequence,
                            accession: psm.ProteinAccession, name: precursorMass.ToString() + " | " + psm.MatchedIonCounts), new DigestionParams(), oneBasedStartResidueInProtein: 1,
                        oneBasedEndResidueInProtein: peptide.BaseSequence.ToCharArray().Length,
                        CleavageSpecificity.Full, peptideDescription: spectrum.OneBasedScanNumber.ToString(), missedCleavages: 1,
                        allModsOneIsNterminus: candidate, numFixedMods: 0);

                    var products = new List<Product>();

                    modifiedPeptide.Fragment(dissociationType: DissociationType.HCD,
                        fragmentationTerminus: FragmentationTerminus.Both, products: products);

                    var match = MetaMorpheusEngine.MatchFragmentIons(new Ms2ScanWithSpecificMass(
                        mzLibScan: spectrum, precursorMonoisotopicPeakMz: spectrum.SelectedIonMZ.Value,
                        precursorCharge: spectrum.SelectedIonChargeStateGuess.Value, fullFilePath: @"D:\08-30-22_bottomup\test.mzML\",
                        new CommonParameters()), products, new CommonParameters());
                    matches.Add(modifiedPeptide, match);
                }

                finalCandidates.Add(matches);
                //finalCandidates.OrderByDescending(x => x.Values.Count);
            }

            return finalCandidates;

        }

        public static double NormalizeArray(double[] array, double currentValue)
        {
            double endOfScale = 1;
            double topOfScale = 100;
            double min = array.Min();
            double max = array.Max();

            var normalized = endOfScale + (currentValue - min) * (topOfScale - endOfScale) / (max - min);

            return normalized;
        }
    }
}

//int aaCounter = 1;
//for (int i = 0; i < baseSequence.ToCharArray().Length; i++)
//{
//    var modsPossibleWithThisResidue = from residue in modDB
//        where residue.Key.ToCharArray()[residue.Key.ToCharArray().Length - 1].Equals(baseSequence.ToCharArray()[i])
//        select residue;

//    if (baseSequence[i].Equals('C'))
//    {
//        if ((modMassCounter + 57.021464).IsSmaller((double)deltaMass, decimalPlaces: 3))
//        {
//            modDictionary.Add(aaCounter, new Modification(_originalId: "Ammonia Loss on C", _monoisotopicMass: -17.026549));
//            modDictionary.Add(aaCounter+1, new Modification(_originalId: "Carbamidomethylation on C", _monoisotopicMass: 57.021464));
//            modMassCounter = modMassCounter + 57.021464 + -17.026549;
//            aaCounter = aaCounter + 2;
//        }
//    }
//    else
//    {
//        foreach (var mod in modsPossibleWithThisResidue)
//        {
//            if ((mod.Value + modMassCounter).IsSmaller((double)deltaMass, 3))
//            {
//                modDictionary.Add(aaCounter, new Modification(_originalId: mod.Key, _monoisotopicMass: mod.Value));
//                modMassCounter = modMassCounter + mod.Value;
//                aaCounter = aaCounter + 1;
//                break;
//            }
//        }
//        aaCounter = aaCounter + 1;
//    }
//}