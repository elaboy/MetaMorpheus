using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Xml;
using Microsoft.Data.Sqlite;

namespace EngineLayer.CombinatorialSearch
{
    public class CSResults : MetaMorpheusEngineResults
    {
        public Dictionary<string, HashSet<Tuple<int, Modification>>> matchedPeptidesDictionary { get; set; }
        public List<Protein> ListOfProteins { get; set; }
        public List<Modification> ModsUsedForSearch { get; set; }

        /// <summary>
        /// Results Constructor.
        /// </summary>
        /// <param name="s"></param>
        /// <param name="matchedPeptidesResults"></param>
        /// <param name="modsUsed"></param>
        public CSResults(MetaMorpheusEngine s, Dictionary<string, HashSet<Tuple<int, Modification>>> matchedPeptidesResults,
            List<Modification> modsUsed, List<Protein> listOfProteins) : base(s)
        {
            matchedPeptidesDictionary = matchedPeptidesResults;
            ModsUsedForSearch = modsUsed;
            ListOfProteins = listOfProteins;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("Proteins trying to expand: " + ListOfProteins.Count);
            sb.AppendLine("Mods types and counts:");
            sb.AppendLine(string.Join(Environment.NewLine, matchedPeptidesDictionary
                    .SelectMany(b => b.Value)
                    .GroupBy(b => b.Item2)
                    .OrderBy(b => -b.Count())
                    .Select(b => "\t" + b.Key.IdWithMotif + "\t" + b.Count())));

            return sb.ToString();
        }
    }
}
