namespace EngineLayer
{
    public class FilteredPsmTSV
    {
        public string FileName;
        public string ScanNumber;
        public string PrecursorScanNumber;
        public string Score;
        public string BaseSeq;
        public string FullSeq;
        public string Mods;
        public string ProteinAccession;
        public string ProteinName;
        public string GeneName;
        public string OrganismName;
        public string StartAndEndResiduesInProtein;
        public string MatchedIonSeries;
        public string MatchedIonCounts;
        public FilteredPsmTSV(string[] psm)
        {
            FileName = psm[0];
            ScanNumber = psm[1];
            PrecursorScanNumber = psm[2];
            Score = psm[3];
            BaseSeq = psm[4];
            FullSeq = psm[5];
            Mods = psm[6];
            ProteinAccession = psm[7];
            ProteinName = psm[8];
            GeneName = psm[9];
            OrganismName = psm[10];
            StartAndEndResiduesInProtein = psm[11];
            MatchedIonSeries = psm[12];
            MatchedIonCounts = psm[13];
        }
    }

}