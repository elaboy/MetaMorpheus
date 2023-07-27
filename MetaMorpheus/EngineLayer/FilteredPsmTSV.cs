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
        public string ModsCount;
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
            ModsCount = psm[7];
            ProteinAccession = psm[8];
            ProteinName = psm[9];
            GeneName = psm[10];
            OrganismName = psm[11];
            StartAndEndResiduesInProtein = psm[12];
            MatchedIonSeries = psm[13];
            MatchedIonCounts = psm[14];
        }
    }

}