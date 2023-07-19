using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using Readers;
using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.IO;
using System.Linq;
using Chemistry;
using EngineLayer;
using MathNet.Numerics;
using TaskLayer;
using Tensorboard;
using ThermoFisher.CommonCore.Data;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class MMGPTMD
    {
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

        public static ImmutableSortedDictionary<double, string> GetModsDictionary()
        {
            var mods = GetModsFromGptmdThing();
            var aaDict = new Mods();

            var modsDictonary = new Dictionary<double, string>();

            foreach (var mod in mods)
            {
                var aminoacid = aaDict.AAsMonoIsotopic.TryGetValue(mod.Target.ToString(), out double residueMass);
                if (residueMass > 0 && modsDictonary.ContainsKey((double)mod.MonoisotopicMass + residueMass) == false)
                {
                    modsDictonary.Add(value: mod.Target.ToString() + " " + mod.OriginalId, key: (double)mod.MonoisotopicMass + residueMass);
                }
            }

            return modsDictonary.ToImmutableSortedDictionary();
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

            var tolerance = new PpmTolerance(30); //Tolerance is better than more or less 

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

                var mods = GetModsDictionary();
                var aaMonoIso = new Mods().AAsMonoIsotopic.Values.ToArray();
                var aaMonoIsoKey = new Mods().AAsMonoIsotopic.Keys.ToArray();
                var ptmMonoIso = mods.Keys.ToArray();
                var ptmMonoIsoResLetter = mods.Values.ToArray();

                foreach (var delta in deltaMz)
                {
                    var possibleSeq =
                        from i in delta
                        where i > 1
                        select i;
                    foreach (var mod in possibleSeq.Select(x => Math.Round(x, 2)).ToList())
                    {
                        //Console.Write("delta: " + mod + " || ");

                        var gotIt = false;
                        for (int i = 0; i < aaMonoIso.Length; i++)
                        {

                            if (mod.Equals(Math.Round(aaMonoIso[i], 2) + 0.01) || mod.Equals(Math.Round(aaMonoIso[i], 2) - 0.01) || mod.Equals(Math.Round(aaMonoIso[i], 2)))
                            {
                                Console.Write(aaMonoIsoKey[i] + " | ");
                                gotIt = true;
                            }
                            if (gotIt)
                            {
                                break;
                            }
                        }

                        for (int k = 0; k < ptmMonoIso.Length; k++)
                        {
                            if (gotIt)
                            {
                                break;
                            }
                            if (mod.Equals(Math.Round(ptmMonoIso[k], 2) + 0.01) || mod.Equals(Math.Round(ptmMonoIso[k], 2) - 0.01) || mod.Equals(Math.Round(ptmMonoIso[k], 2)))
                            {
                                Console.Write(ptmMonoIsoResLetter[k] + " | ");
                                gotIt = true;
                                break;
                            }
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

            //Use list instead of Scans GEScansList()
            for (int i = 0; i < psms.Count(); i++)
            {
                psms[i].PrecursorScanNumber = dataScans[i].OneBasedPrecursorScanNumber.ToString();
                psms[i].ScanNumber = dataScans[i].OneBasedScanNumber.ToString();
            }
            WriteUpdatedFilteredPsmsToTSV(psms, psmFilePath);
        }

        public static void WriteUpdatedFilteredPsmsToTSV(List<FilteredPsmTSV> filteredPsms, string path)
        {
            using (var writer = new StreamWriter(path))
            {
                //makes this an enum?
                string header = "File Name\tScan Number\tPrecursor Scan Number\tScore\tBase Sequence\t" +
                                "Full Sequence\tMods\tProtein Accession\t " +
                                "Protein Name\tGene Name\tOrganism Name\t" +
                                "Start and End Residues in Protein\t" +
                                "Matched Ion Series\tMatched Ion Counts";

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
                        psm.ProteinAccession,
                        psm.ProteinName,
                        psm.GeneName,
                        psm.OrganismName,
                        psm.StartAndEndResiduesInProtein,
                        psm.MatchedIonSeries,
                        psm.MatchedIonCounts
                    };
                    writer.WriteLine(string.Join('\t', row));
                }
            }
        }

        public static void CreateMzML(MsDataScan[] scans, SourceFile sourceFile, string path)
        {
            SourceFile sourceFile1 = new("no nativeID format", "mzML format",
                null, null, filePath: @"K:\08-30-22_bottomup\fractionated_search\Task3-SearchTask\file_example.mzML", null);

            MsDataFile genericFile = new GenericMsDataFile(scans: scans, sourceFile: sourceFile);

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(genericFile, path, true);
        }

        public static (MsDataScan[], MsDataFile) ExtractScansAndSourceFile(List<FilteredPsmTSV> psms, List<string> filePaths)
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
                        var ms1 = file.Value.GetOneBasedScan((int)ms2.OneBasedPrecursorScanNumber);
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
                        MsDataScan scanMS1 = new MsDataScan(massSpectrum: spectrumMS1, oneBasedScanNumber: counter, msnOrder: ms1.MsnOrder, isCentroid: ms1.IsCentroid, polarity: ms1.Polarity,
                            retentionTime: retentionTime, scanWindowRange: ms1.ScanWindowRange, scanFilter: ms1.ScanFilter, mzAnalyzer: ms1.MzAnalyzer,
                            totalIonCurrent: ms1.TotalIonCurrent, injectionTime: injectionTime, noiseData: ms1.NoiseData,
                            nativeId: "controllerType=0 controllerNumber=1 scan=" + counter.ToString(), selectedIonMz: ms1.SelectedIonMZ,
                            selectedIonChargeStateGuess: ms1.SelectedIonChargeStateGuess, selectedIonIntensity: ms1.SelectedIonIntensity,
                            isolationMZ: ms1.IsolationMz, isolationWidth: ms1.IsolationWidth, dissociationType: ms1.DissociationType,
                            hcdEnergy: ms1.HcdEnergy);

                        counter++;
                        retentionTime++;
                        injectionTime++;

                        MsDataScan scanMS2 = new MsDataScan(massSpectrum: spectrumMS2, oneBasedScanNumber: counter, msnOrder: ms2.MsnOrder, isCentroid: ms2.IsCentroid, polarity: ms2.Polarity,
                        retentionTime: retentionTime, scanWindowRange: ms2.ScanWindowRange, scanFilter: ms2.ScanFilter, mzAnalyzer: ms2.MzAnalyzer,
                        totalIonCurrent: ms2.TotalIonCurrent, injectionTime: injectionTime, noiseData: ms2.NoiseData,
                        nativeId: "controllerType=0 controllerNumber=1 scan=" + counter.ToString(), selectedIonMz: ms2.SelectedIonMZ,
                        selectedIonChargeStateGuess: ms2.SelectedIonChargeStateGuess, selectedIonIntensity: ms2.SelectedIonIntensity,
                        isolationMZ: ms2.IsolationMz, isolationWidth: ms2.IsolationWidth, dissociationType: ms2.DissociationType, oneBasedPrecursorScanNumber: precursorNumber,
                        hcdEnergy: ms2.HcdEnergy);


                        counter++;
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

            MsDataScan[] scansArray = scanList.Select(x => x).ToArray();

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
                where psm.QValue <= 0.00001 && psm.PEP <= 0.00001 && psm.DecoyContamTarget.Equals("T") &&
                      PsmFromTsv.ParseModifications(psm.FullSequence).Count() >= 2
                select psm;

            return filteredPsms;
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
                        String.Join('-',psm.FileNameWithoutExtension.Split('-').SkipLast(1)),
                        psm.Ms2ScanNumber.ToString(),
                        psm.PrecursorScanNum.ToString(),
                        psm.RetentionTime.ToString(),
                        string.Join(' ',psm.MatchedIons.Select(x => x.Annotation).ToArray()),
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
                                "Full Sequence\tMods\tProtein Accession\t " +
                                "Protein Name\tGene Name\tOrganism Name\t" +
                                "Start and End Residues in Protein\t" +
                                "Matched Ion Series\tMatched Ion Counts";

                writer.WriteLine(header);
                foreach (var psm in filteredPsms)
                {
                    string[] row = new[]
                    {
                        String.Join('-',psm.FileNameWithoutExtension.Split('-').SkipLast(1)),
                        psm.Ms2ScanNumber.ToString(),
                        psm.PrecursorScanNum.ToString(),
                        psm.Score.ToString(),
                        psm.BaseSeq,
                        psm.FullSequence,
                        PsmFromTsv.ParseModifications(psm.FullSequence).Count().ToString(),
                        psm.ProteinAccession,
                        psm.ProteinName,
                        psm.GeneName,
                        psm.OrganismName,
                        psm.StartAndEndResiduesInProtein,
                        string.Join(' ',psm.MatchedIons.Select(x => x.Annotation).ToArray()),
                        psm.MatchedIons.Count().ToString()
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
}
