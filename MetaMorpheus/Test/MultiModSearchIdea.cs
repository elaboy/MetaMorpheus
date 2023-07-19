using EngineLayer;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Serialization;
using Newtonsoft.Json;
using Proteomics;
using TaskLayer;
using System.Threading.Tasks;

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
}
