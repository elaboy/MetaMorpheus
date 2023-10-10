using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.CombinatorialSearch
{
    public class SwappablePeptideWithSetModifications : PeptideWithSetModifications
    {
        
        public SwappablePeptideWithSetModifications(PeptideWithSetModifications peptide,
            Dictionary<int, Modification> allModsOneIsNterminus) 
            : base(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, "",
                peptide.MissedCleavages, allModsOneIsNterminus, peptide.NumFixedMods, peptide.BaseSequence,
                peptide.PairedTargetDecoyHash)
        {

        }

        public void SwapDict(Dictionary<int, Modification> dict)
        {
            this._allModsOneIsNterminus = dict;
            this._hasChemicalFormulas = null;
            this._monoisotopicMass = null;
            this._mostAbundantMonoisotopicMass = null;
        }
    }
}
