using System.Collections.Generic;
using Proteomics;
using UsefulProteomicsDatabases;
namespace EngineLayer
{
    public class Mods
    {
        public IEnumerable<Modification> ModsKnown { get; set; }
        public int DataScanNumber { get; set; }

        public Mods(string unimodPath)
        {
            ModsKnown = Loaders.LoadUnimod(unimodPath);
        }
    }

}
