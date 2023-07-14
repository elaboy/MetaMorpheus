using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;

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
                Frac1, Frac1_2, Frac2, Frac3, Frac4, Frac5, Frac6, Frac7, Frac8
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
                null, null, filePath: @"D:\08-30-22_bottomup\test.mzML", null));// dataFile.GetSourceFile());
            // todo: update PSMTSV to reflect new mzML file, maybe modify after mzML creation or before?? Before would imply carrying counters maybe as an array or list?


            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msDataFile, @"D:\08-30-22_bottomup\test.mzML", true);

            MMGPTMD.UpdateTheFilteredPsmFile(msDataFile, @"D:\08-30-22_bottomup\example.psmtsv");

            Console.WriteLine("Done");

            Assert.Pass();

        }


    }



}
