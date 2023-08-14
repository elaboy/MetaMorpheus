using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Nett.TomlObjectFactory;

namespace TaskLayer
{
    public class MultipleSearchResults
    {
        public string BaseSequence { get; set; }
        public double SequenceCoverage { get; set; }
        public int IonMatchedCount { get; set; }
        public string Modifications { get; set; }
        public string FullSequence { get; set; }
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

        public static List<DataTable> GetDataTables(List<IGrouping<string, MultipleSearchResults>> results)
        {
            List<DataTable> proteinGroupsTables = new();
            Parallel.ForEach(results, result =>
            {
                var table = new DataTable();
                table.TableName = result.Key;
                foreach (var feature in typeof(MultipleSearchResults).GetProperties())
                {
                    table.Columns.Add(new DataColumn(feature.Name));
                }
                foreach (var peptide in result)
                {
                    DataRow row = table.NewRow();
                    row[0] = peptide.BaseSequence;
                    row[1] = peptide.SequenceCoverage;
                    row[2] = peptide.IonMatchedCount;
                    row[3] = peptide.Modifications;
                    row[4] = peptide.FullSequence;
                    row[5] = peptide.AccessionNumber;
                    row[6] = peptide.PeptideLength;
                    row[7] = peptide.MonoisotopicMass;
                    row[8] = peptide.MostAbundantMonoisotopicMass;
                    row[9] = peptide.IsDecoy;
                    row[10] = String.Join(", ", peptide.MatchedIons);
                    row[11] = String.Join(", ", peptide.MatchedIonCharge);
                    row[12] = String.Join(", ", peptide.TheoricalMz);
                    row[13] = String.Join(", ", peptide.MatchedMz);
                    row[14] = String.Join(", ", peptide.MassErrorPpm);
                    row[15] = String.Join(", ", peptide.MassErrorDa);

                    table.Rows.Add(row);

                }
                proteinGroupsTables.Add(table);
            });
            return proteinGroupsTables;
        }
    }


}
