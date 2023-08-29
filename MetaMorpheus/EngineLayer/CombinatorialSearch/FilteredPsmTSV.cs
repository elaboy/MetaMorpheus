namespace EngineLayer.CombinatorialSearch
{
    public class FilteredPsmTSV
    {
        public string FileName { get; set; }
        public string ScanNumber { get; set; }
        public string PrecursorScanNumber { get; set; }
        public string Score { get; set; }
        public string BaseSeq { get; set; }
        public string FullSeq { get; set; }
        /// <summary>
        /// refers to the variable mods, not fixed mods
        /// </summary>
        public string Mods { get; set; }
        public string ModsCount { get; set; }
        public string ProteinAccession { get; set; }
        public string ProteinName { get; set; }
        public string GeneName { get; set; }
        public string OrganismName { get; set; }
        public string StartAndEndResiduesInProtein { get; set; }
        public string MatchedIonSeries { get; set; }
        public string MatchedIonCounts { get; set; }
        public string PrecursorMass { get; set; }
        public string Charge { get; set; }
        public FilteredPsmTSV(string[] psm)
        {
            FileName = psm[0];
            ScanNumber = psm[1];
            PrecursorScanNumber = psm[2];
            Score = psm[3];
            BaseSeq = psm[4];
            FullSeq = psm[5];
            Mods = psm[6];
            ModsCount = psm[7];
            ProteinAccession = psm[8];
            ProteinName = psm[9];
            GeneName = psm[10];
            OrganismName = psm[11];
            StartAndEndResiduesInProtein = psm[12];
            MatchedIonSeries = psm[13];
            MatchedIonCounts = psm[14];
            PrecursorMass = psm[15];
            Charge = psm[16];
        }

        public FilteredPsmTSV()
        {
        }
    }

}