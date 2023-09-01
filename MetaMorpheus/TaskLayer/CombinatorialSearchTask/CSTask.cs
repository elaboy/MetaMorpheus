using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.CombinatorialSearch;
using MassSpectrometry;
using MzLibUtil;
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
            
            // Classic Search 

            LoadModifications(taskId, out var variableModifications, out var fixedModifications,
                out var localizableModificationTypes);

            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, true, DecoyType.Reverse,
                localizableModificationTypes, CommonParameters);

            List<PeptideSpectralMatch> allPsms = new List<PeptideSpectralMatch>();

            //Engine instance to get mods mass array for the MassDiffAcceptor, will be later running with .Run()

            var engine = new CSEngine(psms, commonBiologicalMods, 3, fixedMods, true,
                new CommonParameters(), dataFile, fileSpecificParameters,
                new List<string>() { "Combinatorial-Search" });
            MyFileManager myFileManager = new MyFileManager(true);

            for (int i = 0; i < currentRawFileList.Count; i++)
            {
                var dataFileName = currentRawFileList[i];

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[i]);

                MassDiffAcceptor searchMode = new DotMassDiffAcceptor("", engine.MassArray.Distinct().AsEnumerable(),
                    combinedParams.PrecursorMassTolerance);

                NewCollection(Path.GetFileName(dataFileName), new List<string>{taskId, "Individual Spectra Files", dataFileName});
                MsDataFile myDataFile = myFileManager.LoadFile(dataFileName, combinedParams);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass =
                    GetMs2Scans(myDataFile, dataFileName, combinedParams).OrderBy(x => x.PrecursorMass).ToArray();
                myFileManager.DoneWithFile(dataFileName);
                PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[arrayOfMs2ScansSortedByMass.Length];

                // search
                new ClassicSearchEngine(allPsmsArray, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications,
                    null, null, null,
                    proteinList, searchMode, combinedParams, this.FileSpecificParameters, null,
                    new List<string>(){taskId, "Individual Spectra Files", dataFileName}, false).Run();

                allPsms.AddRange(allPsmsArray.Where(x => x != null));
                FinishedDataFile(dataFileName, new List<string>(){taskId, "Individual Spectra Files", dataFileName});
            }

            allPsms = allPsms.OrderByDescending(b => b.Score)
                .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue)
                .GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass))
                .Select(b => b.First()).ToList();

            // CS stuff 
            var csResults = (CSResults)engine.Run();


            
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
