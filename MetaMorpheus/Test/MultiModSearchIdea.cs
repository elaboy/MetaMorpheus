using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Readers;
using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.IO;
using System.Linq;
using System.Threading;
using System.Windows.Automation.Peers;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class MultiModSearchIdea
    {
        #region Files Path

        private const string Frac1 =
            @"D:\\08-30-22_bottomup\\fractionated\08-31-22_fractionated_human_Tryp_40ug_F1.raw";

        private const string Frac1_2 =
            @"D:\\08-30-22_bottomup\\fractionated\\08-31-22_fractionated_human_Tryp_40ug_F1_2.raw";

        private const string Frac2 =
            @"D:\\08-30-22_bottomup\\fractionated\\08-31-22_fractionated_human_Tryp_40ug_F2.raw";

        private const string Frac3 =
            @"D:\\08-30-22_bottomup\\fractionated\\08-31-22_fractionated_human_Tryp_40ug_F3.raw";

        private const string Frac4 =
            @"D:\\08-30-22_bottomup\\fractionated\\08-31-22_fractionated_human_Tryp_40ug_F4.raw";

        private const string Frac5 =
            @"D:\\08-30-22_bottomup\\fractionated\\08-31-22_fractionated_human_Tryp_40ug_F5.raw";

        private const string Frac6 =
            @"D:\\08-30-22_bottomup\\fractionated\\08-31-22_fractionated_human_Tryp_40ug_F6.raw";

        private const string Frac7 =
            @"D:\\08-30-22_bottomup\\fractionated\\08-31-22_fractionated_human_Tryp_40ug_F7.raw";

        private const string Frac8 =
            @"D:\\08-30-22_bottomup\\fractionated\\08-31-22_fractionated_human_Tryp_40ug_F8.raw";

        private const string AllPeptidesPsm1 =
            @"D:\08-30-22_bottomup\fractionated_search\Task3-SearchTask\AllPeptides.psmtsv";
        #endregion

        [Test]
        public void Run()
        {
            List<string> filePaths = new List<string>()
            {
                Frac1, Frac1_2, Frac3, Frac4, Frac5, Frac6, Frac7, Frac8
            };

            var psms = PsmTsvReader.ReadTsv(AllPeptidesPsm1, out List<string> warnings);

            var filteredPsms = FilterPsm(psms).ToList();

            WriteFilteredPsmsToTSV(filteredPsms, @"D:\08-30-22_bottomup\fractionated_search\Task3-SearchTask\example.psmtsv");

            var filteredFile =
                ReadFilteredPsmTSV(@"D:\08-30-22_bottomup\fractionated_search\Task3-SearchTask\example.psmtsv");

            WriteFastaDBFromFilteredPsm(filteredFile,
                @"D:\08-30-22_bottomup\fractionated_search\Task3-SearchTask\database_example.fasta");

            var scans = ExtractScans(filteredFile, filePaths);

            CreateMzML(scans, @"D:\08-30-22_bottomup\fractionated_search\Task3-SearchTask\file_example.mzML");

            Console.WriteLine("Done");

        }

        public static void CreateMzML(List<MsDataScan> scans, string path)
        {
            SourceFile sourceFile = new SourceFile("no nativeID format", "mzML format",
                null, null, null);

            //var file = new GenericMsDataFile(scans.ToArray(), sourceFile);
            //File.Create(path);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(scans.ToArray(), sourceFile), path, false);
        }

        public List<MsDataScan> ExtractScans(List<FilteredPsmTSV> psms, List<string> filePaths)
        {
            List<MsDataFile> loadedFiles = new();

            foreach (var file in filePaths)
            {
                int worker = 0;
                int io = 0;
                ThreadPool.GetAvailableThreads(workerThreads: out worker, completionPortThreads: out io);

                loadedFiles.Add(new ThermoRawFileReader(file).LoadAllStaticData(maxThreads: worker));
            }

            List<string> fileName = new();

            foreach (var name in loadedFiles)
            {
                fileName.Add(name.FilePath);
            }

            List<Tuple<string, MsDataFile>> tupleList = new();

            for (int i = 0; i < loadedFiles.Count(); i++)
            {
                tupleList.Add(new Tuple<string, MsDataFile>(item1: fileName[i], item2: loadedFiles[i]));
            }

            var dict = tupleList.ToImmutableDictionary(x => x.Item1, x => x.Item2);

            List<MsDataScan> scans = new();

            foreach (var psm in psms)
            {
                foreach (var item in dict)
                {
                    if (item.Value.SourceFile.FileName.Contains(psm.FileName))
                    {
                        foreach (var scan in item.Value.Scans)
                        {
                            if (scan.OneBasedScanNumber.Equals(int.Parse(psm.ScanNumber)) && scan.MsnOrder == 2)
                            {
                                scans.Add(scan);
                                break;
                            }
                        }

                        break;
                    }

                    //break;

                    //if (item.Value.SourceFile.FileName.Contains(psm.FileName) && item.Value.NumSpectra.Equals(int.Parse(psm.ScanNumber)))
                    //{
                    //    scans.Add(item.Value.GetOneBasedScan(int.Parse(psm.ScanNumber)));
                    //    break;
                    //}
                }
            }

            //ImmutableDictionary<string, MsDataFile> files =
            //    ImmutableDictionary.CreateRange(
            //        new KeyValuePair<string, MsDataFile>[]
            //        {
            //            KeyValuePair.Create("", new ThermoRawFileReader(filePath[0]).LoadAllStaticData())
            //        }
            //    );

            //foreach (var path in filePath)
            //{
            //    files.Add(new ThermoRawFileReader(path));
            //}

            //foreach (var psm in psms)
            //{
            //    scans.Add(files.GetOneBasedScan(int.Parse(psm.ScanNumber)));
            //}

            return scans;
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


        public static void WriteFilteredPsmsToTSV(List<PsmFromTsv> filteredPsms, string path)
        {
            using (var writer = new StreamWriter(path))
            {
                //makes this an enum?
                string header = "File Name\tScan Number\tScore\tBase Sequence\t" +
                                "Full Sequence\tMods\tProtein Accession\t " +
                                "Protein Name\tGene Name\tOrganism Name\t" +
                                "Start and End Residues in Protein\t" +
                                "Matched Ion Series\tMatched Ion Counts";
                writer.WriteLine(header);
                foreach (var psm in filteredPsms)
                {
                    string[] row = new[]
                    {
                        String.Join('_',psm.FileNameWithoutExtension.Split('_').SkipLast(1)),
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

        public List<FilteredPsmTSV> ReadFilteredPsmTSV(string path)
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

    public class FilteredPsmTSV
    {
        public string FileName;
        public string ScanNumber;
        public string PrecursorScanNumber;
        public string Score;
        public string BaseSeq;
        public string FullSeq;
        public string Mods;
        public string ProteinAccession;
        public string ProteinName;
        public string GeneName;
        public string OrganismName;
        public string StartAndEndResiduesInProtein;
        public string MatchedIonSeries;
        public string MatchedIonCounts;

        public FilteredPsmTSV(string[] psm)
        {
            FileName = psm[0];
            ScanNumber = psm[1];
            PrecursorScanNumber = psm[2];
            Score = psm[3];
            BaseSeq = psm[4];
            FullSeq = psm[5];
            Mods = psm[6];
            ProteinAccession = psm[7];
            ProteinName = psm[8];
            GeneName = psm[9];
            OrganismName = psm[10];
            StartAndEndResiduesInProtein = psm[11];
            MatchedIonSeries = psm[12];
            MatchedIonCounts = psm[13];
        }
    }

}
