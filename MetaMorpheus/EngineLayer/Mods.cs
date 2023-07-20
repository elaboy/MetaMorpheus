using System.Collections.Generic;
using System.Collections.Immutable;
using Proteomics;
using UsefulProteomicsDatabases;
namespace EngineLayer
{
    public class Mods
    {
        public IEnumerable<Modification> ModsKnown { get; set; }
        public int DataScanNumber { get; set; }

        // http://www.matrixscience.com/help/aa_help.html
        public Dictionary<string, double> AAsMonoIsotopic = new Dictionary<string, double>(
            new KeyValuePair<string, double>[]
            {
                KeyValuePair.Create("A", 71.037114),
                KeyValuePair.Create("R", 156.101111),
                KeyValuePair.Create("N", 114.042927),
                KeyValuePair.Create("D", 115.026943),
                KeyValuePair.Create("C", 103.009185),
                KeyValuePair.Create("E", 129.042593),
                KeyValuePair.Create("Q", 128.058578),
                KeyValuePair.Create("G", 57.021464),
                KeyValuePair.Create("H", 137.058912),
                KeyValuePair.Create("I", 113.084064),
                KeyValuePair.Create("L", 113.084064),
                KeyValuePair.Create("K", 128.094963),
                KeyValuePair.Create("M", 131.040485),
                KeyValuePair.Create("F", 147.068414),
                KeyValuePair.Create("P", 97.052764),
                KeyValuePair.Create("S", 87.032028),
                KeyValuePair.Create("T", 101.047679),
                KeyValuePair.Create("U", 150.95363),
                KeyValuePair.Create("W", 186.079313),
                KeyValuePair.Create("Y", 163.06332),
                KeyValuePair.Create("V", 99.068414)
            });

        public Mods()
        {

        }
        public Mods(string unimodPath)
        {
            ModsKnown = Loaders.LoadUnimod(unimodPath);
        }


    }

}
