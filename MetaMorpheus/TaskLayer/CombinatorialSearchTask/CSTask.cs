using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using EngineLayer.CombinatorialSearch;
using Proteomics;
using UsefulProteomicsDatabases;

namespace TaskLayer.CombinatorialSearchTask
{
    public class CSTask : MetaMorpheusTask
    {
        public CommonParameters CommonParameters { get; set; }
        public CSTask() : base(MyTask.CombinatorialSearch)
        {

        }

        protected override MyTaskResults RunSpecific(string OutputFolder,
            List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId,
            FileSpecificParameters[] fileSettingsList)
        {
            // setup 

            var mods =
                Loaders.LoadUnimod(
                    @"C:\Users\elabo\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml").ToList();

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


            // search 

            // CS stuff 
            var engine = new CSEngine(psms, commonBiologicalMods, 4, fixedMods, true,
                new CommonParameters(), @"D:\08-30-22_bottomup\test.mzML");

            // output results


            throw new NotImplementedException();
        }
    }
}
