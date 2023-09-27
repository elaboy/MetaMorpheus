using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.CombinatorialSearch;
using MassSpectrometry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Data;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using EngineLayer.FdrAnalysis;
using Microsoft.ML.Transforms;
using Proteomics.Fragmentation;
using Readers;
using UsefulProteomicsDatabases;

namespace TaskLayer.CombinatorialSearchTask
{
    public class CSTask : MetaMorpheusTask
    {
        private Dictionary<string, int> XMLThing { get; set; }

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

            //var psms = CSEngine.ReadFilteredPsmTSVShort(@"D:\08-30-22_bottomup\example.psmtsv"); //not being used rn

            //var dataFile = Readers.MsDataFileReader.GetDataFile(@"D:\08-30-22_bottomup\test.mzML");
            List<(string, CommonParameters)> fileSpecificParameters = new();

            foreach (var filePath in currentRawFileList)
            {
                var dataFile = MsDataFileReader.GetDataFile(filePath);
                fileSpecificParameters.Add((dataFile.GetSourceFile().FileName,
                    CommonParameters));
            } //todo load this inside the forLoop


            MyTaskResults = new(this) { NewDatabases = new List<DbForTask>() };

            // Classic Search 

            LoadModifications(taskId, out var variableModifications, out var fixedModifications,
                out var localizableModificationTypes);

            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, true, DecoyType.Reverse,
                localizableModificationTypes, CommonParameters);

            List<PeptideSpectralMatch> allPsms = new List<PeptideSpectralMatch>();

            //Engine instance to get mods mass array for the MassDiffAcceptor, will be later running with .Run()
            List<List<Modification>> listOfModCombinations = new();

            CSExtension.CombinationBuilder(commonBiologicalMods, ref listOfModCombinations, 2, true);

            CSExtension.BuildCombinationWithAddedMass(listOfModCombinations, out var combinationsWithAddedMass);


            //var engine = new CSEngine(psms, commonBiologicalMods, 3, fixedMods, true,
            //    new CommonParameters(), dataFile, fileSpecificParameters,
            //    new List<string>() { "Combinatorial-Search" });

            MyFileManager myFileManager = new MyFileManager(true);
            for (int i = 0; i < currentRawFileList.Count; i++)
            //Parallel.For(0, currentRawFileList.Count, i =>
            {
                var dataFileName = currentRawFileList[i];

                CommonParameters combinedParams =
                    SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[i]);

                MassDiffAcceptor searchMode = new DotMassDiffAcceptor("",
                    commonBiologicalMods.Select(x => x.MonoisotopicMass.Value).Distinct(),
                    combinedParams.PrecursorMassTolerance);

                NewCollection(Path.GetFileName(dataFileName),
                    new List<string> { taskId, "Individual Spectra Files", dataFileName });
                MsDataFile myDataFile = myFileManager.LoadFile(dataFileName, combinedParams);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass =
                    GetMs2Scans(myDataFile, dataFileName, combinedParams)
                        .OrderBy(x => x.PrecursorMass).ToArray();

                myFileManager.DoneWithFile(dataFileName);

                PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[arrayOfMs2ScansSortedByMass.Length];

                // search
                new ClassicSearchEngine(allPsmsArray, arrayOfMs2ScansSortedByMass, variableModifications,
                    new List<Modification>(),
                    null, null, null,
                    proteinList, searchMode, combinedParams, this.FileSpecificParameters, null,
                    new List<string>() { taskId, "Individual Spectra Files", dataFileName }, false).Run();

                allPsms.AddRange(allPsmsArray.Where(x => x != null));
                FinishedDataFile(dataFileName,
                    new List<string>() { taskId, "Individual Spectra Files", dataFileName });
            }


            allPsms = allPsms
                .OrderByDescending(b => b.Score)
                .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue)
                .GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass))
                .Select(b => b.First()).ToList();

            MassDiffAcceptor searchModeTemp = new DotMassDiffAcceptor("", commonBiologicalMods
                    .Select(x => x.MonoisotopicMass.Value), CommonParameters.PrecursorMassTolerance);

            new FdrAnalysisEngine(allPsms, searchModeTemp.NumNotches, CommonParameters,
                this.FileSpecificParameters, new List<string>() { taskId });

            WritePsmsToTsv(allPsms, Path.Combine(OutputFolder, " CS_Candidates.psmtsv"), new Dictionary<string, int>());
            FinishedWritingFile(Path.Combine(OutputFolder, " CS_Candidates.psmtsv"), new List<string> { taskId });

            // CS stuff 


            var engine = new CSEngine(allPsms
                    .Where(x => x.IsDecoy == false && x.BaseSequence is not null), proteinList, 
                listOfModCombinations, combinationsWithAddedMass, new List<Modification>(), //changed the fixed mods to treat as variable in test case
                new CommonParameters(), fileSpecificParameters,
                new List<string>() { "Combinatorial-Search" });

            var csResults = (CSResults)engine.Run();

            ProteinDbWriter.WriteXmlDatabase(csResults.matchedPeptidesDictionary,
                proteinList, @"D:\TestingCSTask\testingTasks.xml");

            FinishedWritingFile(@"D:\TestingCSTask", new List<string>() { taskId });
            MyTaskResults.NewDatabases.Add(new DbForTask(@"D:\08-30-22_bottomup\test.mzML", false));

            return MyTaskResults;
        }

    }
}
