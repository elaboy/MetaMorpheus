using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using EngineLayer.CombinatorialSearch;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using ScottPlot.Statistics;
using UsefulProteomicsDatabases;

namespace TaskLayer.CombinatorialSearchTask
{
    public class CSTask : MetaMorpheusTask
    {
        public CSTask(CommonParameters commonParameters) : base(MyTask.CombinatorialSearch)
        {
            CommonParameters = commonParameters;
        }

        protected override MyTaskResults RunSpecific(string OutputFolder,
            List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId,
            FileSpecificParameters[] fileSettingsList)
        {
            // setup 

            var mods =
                Loaders.LoadUnimod(
                    @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml").ToList();

            var fixedMods = new List<Modification>();
            fixedMods.Add(mods.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")));

            List<Modification> commonBiologicalMods = GlobalVariables.AllModsKnown.OfType<Modification>()
                .Where(mod =>
                    mod.ModificationType.Contains("Common Biological") ||
                    mod.ModificationType.Contains("Common Artifact") ||
                    mod.ModificationType.Contains("Less Common") ||
                    mod.ModificationType.Contains("Metals") ||
                    mod.ModificationType.Contains("UniProt")).ToList();

            var psms = CSEngine.ReadFilteredPsmTSVShort(@"D:\08-30-22_bottomup\example.psmtsv");

            var dataFile = Readers.MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML");
            List<(string, CommonParameters)> fileSpecificParameters = new();

            fileSpecificParameters.Add((dataFile.GetSourceFile().FileName,
                CommonParameters));

            MyTaskResults = new(this) { NewDatabases = new List<DbForTask>() };
            // search 

            // CS stuff 
            var csResults = (CSResults)new CSEngine(psms, commonBiologicalMods, 3, fixedMods, true,
                new CommonParameters(), dataFile, fileSpecificParameters, new List<string>(){ "Combinatorial-Search" })
                .Run();


            
            // output results
            //var writtenModsIntoDB = ProteinDbWriter.WriteXmlDatabase()
            var xmlDatabase = ProteinDbWriter.WriteXmlDatabase(csResults.matchedPeptidesDictionary,
                csResults.ListOfProteins, @"D:\TestingCSTask\testingTasks.xml");

            FinishedWritingFile(@"D:\TestingCSTask", new List<string>() {taskId});
            MyTaskResults.NewDatabases.Add(new DbForTask(@"D:\08-30-22_bottomup\test.mzML", false));

            return MyTaskResults;
        }
    }
}
