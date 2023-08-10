using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using EngineLayer;
using MassSpectrometry;
using Microsoft.Win32;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Readers;
using TaskLayer;
using Test;
using UsefulProteomicsDatabases;
using static Test.TestModComboSearch;

namespace MultiModSearchViz
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void LoadPsm(object sender, RoutedEventArgs e)
        {
            OpenFileDialog psmFile = new OpenFileDialog();
            psmFile.ShowDialog();
            psmFile.FileName = @"Asdad";

            var psmList = MultiModSearch.ReadFilteredPsmTSVShort(psmFile.FileName);

        }

        private void LoadDataFile(object sender, RoutedEventArgs e)
        {
            OpenFileDialog dataFilePath = new OpenFileDialog();
            dataFilePath.ShowDialog();
            dataFilePath.FileName = @"asdada";

            var dataFile =MsDataFileReader.GetDataFile(dataFilePath.FileName).GetAllScansList()
                .Where(x => x.MsnOrder == 2).ToArray();
        }

        private void RunButton(object sender, RoutedEventArgs e)
        {
            List<FilteredPsmTSV> psmList = MultiModSearch.ReadFilteredPsmTSVShort($@"{psmPath.Text}");

            MsDataFile dataFile = MsDataFileReader.GetDataFile($@"{mzMLPath.Text}");

            List<Modification> commonBiologycalMods = GlobalVariables.AllModsKnown.OfType<Modification>()
                .Where(mod => mod.ModificationType.Contains("Common Biological")).ToList();

            List<Modification> modsFromUnimod =
                Loaders.LoadUnimod(
                        @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                    .ToList();

            List<Modification> fixedMods = new List<Modification>();
            fixedMods.Add(modsFromUnimod.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")));

            MultipleSearchEngine engine = new MultipleSearchEngine(commonBiologycalMods, int.Parse(maxNumOfMods.Text), true);

            List<KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>> tempMatches = new();


            Parallel.ForEach(psmList, psm =>
            {
                IEnumerable<PeptideWithSetModifications> peptideAsProtein = MultipleSearchEngine.GetPeptideAsProtein(psm, fixedMods);
                List<KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>>> matches = 
                    engine.GetPeptideFragmentIonsMatches(psm, dataFile, fixedMods);

                tempMatches.AddRange(matches);
            });

            var results = tempMatches.OrderByDescending(x => x.Value.Count)
                .GroupBy(x => x.Key.BaseSequence);


            List<DataTable> proteinGroupsTables = new();

            Parallel.ForEach(results, result =>
            {
                var table = new DataTable();
                foreach (var feature in typeof(MultiModSearchResults).GetProperties())
                {
                    table.Columns.Add(new DataColumn(feature.Name));
                }
                foreach (var peptide in result)
                {
                    var individualPeptide = new MultiModSearchResults()
                    {
                        AccessionNumber = peptide.Key.Protein.Accession,
                        BaseSequece = peptide.Key.Protein.BaseSequence,
                        IsDecoy = peptide.Key.Protein.IsDecoy,
                        MassErrorDa = peptide.Value.Select(x => x.MassErrorDa).ToArray(),
                        MassErrorPpm = peptide.Value.Select(x => x.MassErrorPpm).ToArray(),
                        MatchedIonCharge = peptide.Value.Select(x => x.Charge).ToArray(),
                        MatchedIons = peptide.Value.Select(x => x.Annotation).ToArray(),
                        MatchedMz = peptide.Value.Select(x => x.Mz).ToArray(),
                        TheoricalMz = peptide.Value.Select(x => x.NeutralTheoreticalProduct.NeutralMass).ToArray(),
                        MonoisotopicMass = peptide.Key.MonoisotopicMass,
                        MostAbundantMonoisotopicMass = peptide.Key.MostAbundantMonoisotopicMass,
                        PeptideLength = peptide.Key.Protein.Length
                    };
                    DataRow row = table.NewRow();
                    row[1] = individualPeptide.BaseSequece;
                    row[2] = individualPeptide.AccessionNumber;
                    row[3] = individualPeptide.PeptideLength;
                    row[3] = individualPeptide.MonoisotopicMass;
                    row[4] = individualPeptide.MostAbundantMonoisotopicMass;
                    row[5] = individualPeptide.IsDecoy;
                    row[6] = String.Join(", ", individualPeptide.MatchedIons);
                    row[7] = String.Join(", ", individualPeptide.MatchedIonCharge);
                    row[8] = String.Join(", ", individualPeptide.TheoricalMz);
                    row[9] = String.Join(", ", individualPeptide.MatchedMz);
                    row[10] = String.Join(", ", individualPeptide.MassErrorPpm);
                    row[11] = String.Join(", ", individualPeptide.MassErrorDa);
                }
                proteinGroupsTables.Add(table);
            });

            proteinGroups.ItemsSource = proteinGroupsTables.AsEnumerable();
        }
    }
}
