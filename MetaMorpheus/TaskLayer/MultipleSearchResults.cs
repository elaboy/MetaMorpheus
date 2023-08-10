using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer
{
    public class MultipleSearchResults
    {
        public string BaseSequece { get; set; }
        public string AccessionNumber { get; set; }
        public int PeptideLength { get; set; }
        public double MonoisotopicMass { get; set; }
        public double MostAbundantMonoisotopicMass { get; set; }
        public bool IsDecoy { get; set; }
        public string[] MatchedIons { get; set; }
        public int[] MatchedIonCharge { get; set; }
        public double[] TheoricalMz { get; set; }
        public double[] MatchedMz { get; set; }
        public double[] MassErrorPpm { get; set; }
        public double[] MassErrorDa { get; set; }
    }
}
