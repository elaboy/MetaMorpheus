using EngineLayer;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Xml.Serialization;
using Newtonsoft.Json;
using Proteomics;
using TaskLayer;
using System.Threading.Tasks;
using iText.IO.Source;
using iText.Kernel.Pdf.Canvas.Parser.ClipperLib;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;
using ScottPlot;
using ScottPlot.Renderable;
using TorchSharp;
using MathNet.Numerics.Random;

namespace Test
{
    [TestFixture]
    public class MultiModSearchIdea
    {
        #region Files Path

        private const string Frac1 =
            @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F1-calib.mzML";

        //private const string Frac1_2 =
        //    @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F2-calib.mzML";

        private const string Frac2 =
            @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F2-calib.mzML";

        private const string Frac3 =
            @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F3-calib.mzML";

        private const string Frac4 =
            @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F4-calib.mzML";

        private const string Frac5 =
            @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F5-calib.mzML";

        private const string Frac6 =
            @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F6-calib.mzML";

        private const string Frac7 =
            @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F7-calib.mzML";

        private const string Frac8 =
            @"D:\08-30-22_bottomup\fractionated_search\Task1-CalibrateTask\08-31-22_fractionated_human_Tryp_40ug_F8-calib.mzML";

        private const string AllPeptidesPsm1 =
            @"D:\08-30-22_bottomup\fractionated_search\Task3-SearchTask\AllPeptides.psmtsv";

        #endregion

        [Test]
        public void TestMzMLWitter()
        {
            var file = Readers.ThermoRawFileReader.LoadAllStaticData(Frac1);

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file, "test.mzML", false);

            Console.WriteLine("yes");
        }

        [Test]
        public void Run()
        {
            //MsDataFile msDataFile = MsDataFileReader.GetDataFile(@"D:\\08-30-22_bottomup\\test.mzML");
            //msDataFile.LoadAllStaticData();
            //var scanList = msDataFile.GetAllScansList();

            List<string> filePaths = new List<string>()
            {
                Frac1, Frac2, Frac3, Frac4, Frac5, Frac6, Frac7, Frac8
            };

            var psms = PsmTsvReader.ReadTsv(AllPeptidesPsm1, out List<string> warnings);

            var filteredPsms = MMGPTMD.FilterPsm(psms).ToList();



            MMGPTMD.WriteFilteredPsmsToTSV(filteredPsms, @"D:\08-30-22_bottomup\example.psmtsv");

            var filteredFile =
                MMGPTMD.ReadFilteredPsmTSVShort(@"D:\08-30-22_bottomup\example.psmtsv");

            MMGPTMD.WriteFastaDBFromFilteredPsm(filteredFile,
                @"D:\08-30-22_bottomup\database_example.fasta");

            var (scans, dataFile) = MMGPTMD.ExtractScansAndSourceFile(filteredFile, filePaths);

            MsDataFile msDataFile = new GenericMsDataFile(scans, new("no nativeID format", "mzML format",
                null, null, filePath: @"D:\08-30-22_bottomup\test.mzML", null)); // dataFile.GetSourceFile());
            // todo: update PSMTSV to reflect new mzML file, maybe modify after mzML creation or before?? Before would imply carrying counters maybe as an array or list?


            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msDataFile, @"D:\08-30-22_bottomup\test.mzML", false);

            MMGPTMD.UpdateTheFilteredPsmFile(msDataFile, @"D:\08-30-22_bottomup\example.psmtsv");

            Console.WriteLine("Done");

            Assert.Pass();

        }

        [Test]
        public void RunTopDown()
        {
            //MsDataFile msDataFile = MsDataFileReader.GetDataFile(@"D:\\08-30-22_bottomup\\test.mzML");
            //msDataFile.LoadAllStaticData();
            //var scanList = msDataFile.GetAllScansList();

            List<string> filePaths = new List<string>()
            {
                @"D:\topDown\02-18-20_jurkat_td_rep2_fract10.raw",
                @"D:\topDown\02-17-20_jurkat_td_rep2_fract2.raw",
                @"D:\topDown\02-17-20_jurkat_td_rep2_fract3.raw",
                @"D:\topDown\02-17-20_jurkat_td_rep2_fract4.raw",
                @"D:\topDown\02-18-20_jurkat_td_rep2_fract5.raw",
                @"D:\topDown\02-18-20_jurkat_td_rep2_fract6.raw",
                @"D:\topDown\02-18-20_jurkat_td_rep2_fract7.raw",
                @"D:\topDown\02-18-20_jurkat_td_rep2_fract8.raw",
                @"D:\topDown\02-18-20_jurkat_td_rep2_fract9.raw"
            };

            var psms = PsmTsvReader.ReadTsv(@"D:\topDown\MOxAndBioMetArtModsGPTMD_Search\Task2-SearchTask\AllProteoforms.psmtsv", out List<string> warnings);

            var filteredPsms = MMGPTMD.FilterPsm(psms).ToList();



            MMGPTMD.WriteFilteredPsmsToTSV(filteredPsms, @"D:\topDown\example.psmtsv");

            var filteredFile =
                MMGPTMD.ReadFilteredPsmTSVShort(@"D:\topDown\example.psmtsv");

            MMGPTMD.WriteFastaDBFromFilteredPsm(filteredFile,
                @"D:\topDown\database_example.fasta");

            var (scans, dataFile) = MMGPTMD.ExtractScansAndSourceFile(filteredFile, filePaths);

            MsDataFile msDataFile = new GenericMsDataFile(scans, new("no nativeID format", "mzML format",
                null, null, filePath: @"D:\topDown\test.mzML", null)); // dataFile.GetSourceFile());
            // todo: update PSMTSV to reflect new mzML file, maybe modify after mzML creation or before?? Before would imply carrying counters maybe as an array or list?


            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msDataFile, @"D:\topDown\test.mzML", false);

            MMGPTMD.UpdateTheFilteredPsmFile(msDataFile, @"D:\topDown\example.psmtsv");

            Console.WriteLine("Done");

            Assert.Pass();

        }

        [Test]
        public void TestMultiModDiscovery()
        {
            MMGPTMD.MultiModDiscovery(@"D:\08-30-22_bottomup\test.mzML");

            Console.WriteLine("Okay");
        }

        [Test]
        public void TestGetModsFromGptmdThing()
        {
            var test = MMGPTMD.GetModsFromGptmdThing(@"Task1-GPTMDTaskconfig.toml");
        }

        [Test]
        public void TestGetModsDictionary()
        {
            var dict = MMGPTMD.GetModsDictionary();

        }
        [Test]
        public void TestGetCandidates()
        {
            var dict = MMGPTMD.GetCandidates();

        }

        [Test]
        public void TestMatchSpectra()
        {
            MMGPTMD.MatchSpectra();
        }

        [Test]
        public void MMGPTMDRun()
        {
            var search = new MMGPTMD(21);
        }

        [Test]
        public void MultiMMGPTMD()
        {
            List<MMGPTMD> searches = new List<MMGPTMD>();

            for (int i = 0; i < 26; i++)
            {
                searches.Add(new MMGPTMD(i));
            }
        }


        //////[Test]
        //////public void MakeModsDB()
        //////{
        //////    var task = Toml.ReadFile<GptmdTask>(@"C:\Users\Edwin\Downloads\Task1-GPTMDTaskconfig.toml",
        //////        MetaMorpheusTask.tomlConfig);

        //////    var mods = GlobalVariables.AllModsKnownDictionary;

        //////    List<Modification> modificationList = new();

        //////    foreach (var (item1, item2) in task.GptmdParameters.ListOfModsGptmd)
        //////    {
        //////        if (mods.TryGetValue(item2, out Modification mod))
        //////        {
        //////            modificationList.Add(mod);
        //////        }
        //////    }


        //////    XmlSerializer temp = new XmlSerializer(typeof(List<Modification>));
        //////    using (StreamWriter writer = new(File.Create(@"C:\Users\Edwin\Downloads\Mods.xml")))
        //////    {
        //////        temp.Serialize(writer, modificationList);
        //////    }
        //////}

    }

    public sealed class TestScalableModSearch
    {
        private const string _filteredPsm = @"D:\08-30-22_bottomup\example.psmtsv";

        [Test]
        public void TestReadFilteredPsmTSVShort()
        {
            var psmList =
                ScalableModSearch.ReadFilteredPsmTSVShort(_filteredPsm);

            List<Dictionary<int, Modification>> testing = new();
            var msDataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML");

            ScalableModSearch.NewGetAllComboMods(msDataFile, psmList);
            //.GetTopScoreAndSavePng(pathToSavePlots:@"D:/08-30-22_bottomup/plotImages/");


        }

        [Test]
        public void TestSetCandidateList()
        {
            var psmList = ScalableModSearch.ReadFilteredPsmTSVShort(_filteredPsm);
            var msDataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML");

            var mods =
                Loaders.LoadUnimod(
                    @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml");

            var search = new ScalableModSearch(msDataFile, psmList);

            Console.WriteLine();
        }

        [Test]
        public void TestMassSpectrumPlot()
        {
            var msDataFile = MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML");
            MsDataScan scan = msDataFile.GetOneBasedScan(2);

            double[] mz = scan.MassSpectrum.XArray;
            double[] intensity = scan.MassSpectrum.YArray;

            var plt = new ScottPlot.Plot(1000, 700);

            for (int i = 0; i < mz.Length; i++)
            {
                var vlines = new ScottPlot.Plottable.VLineVector();
                vlines.Xs = new[] { mz[i] };
                vlines.Max = intensity[i];
                vlines.Color = Color.Black;
                //vlines.PositionLabel = true;
                //vlines.PositionLabelBackground = vlines.Color;
                plt.Add(vlines);
            }
            plt.SetAxisLimitsY(0, intensity.Max());
            plt.SaveFig("msScanTest.png");

        }

        [Test]
        public void TestVerticalLines()
        {
            var plt = new ScottPlot.Plot(600, 400);

            Random rand = new Random(0);
            double[] xs = DataGen.Random(rand, 50);
            double[] ys = DataGen.Random(rand, 50);

            //var scatter = plt.AddScatterPoints(xs, ys, Color.Blue, 10);

            var vlines = new ScottPlot.Plottable.VLineVector();
            vlines.Xs = new double[] { xs[1], xs[12], xs[35] };
            vlines.Color = Color.Red;
            vlines.PositionLabel = true;
            vlines.PositionLabelBackground = vlines.Color;

            var hlines = new ScottPlot.Plottable.HLineVector();
            hlines.Ys = new double[] { ys[1], ys[12], ys[35] };
            hlines.Color = Color.DarkCyan;
            hlines.PositionLabel = true;
            hlines.PositionLabelBackground = hlines.Color;
            hlines.DragEnabled = true;

            //plt.Add(scatter);
            plt.Add(vlines);
            plt.Add(hlines);

            plt.SaveFig("axisLine_Vector.png");
        }
    }
}
