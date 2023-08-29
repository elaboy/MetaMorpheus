using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;

namespace EngineLayer.CombinatorialSearch
{
    public class CSResults : MetaMorpheusEngineResults
    {
        private Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> matchedPeptidesDictionary { get; set; }
        private List<Modification> ModsUsedForSearch { get; set; }

        /// <summary>
        /// Results Constructor.
        /// </summary>
        /// <param name="s"></param>
        /// <param name="matchedPeptidesResults"></param>
        /// <param name="modsUsed"></param>
        public CSResults(MetaMorpheusEngine s, Dictionary<PeptideWithSetModifications,
            List<MatchedFragmentIon>> matchedPeptidesResults, List<Modification> modsUsed) : base(s)
        {
            matchedPeptidesDictionary = matchedPeptidesResults;
            ModsUsedForSearch = modsUsed;

        }
    }
}
