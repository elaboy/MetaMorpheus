using EngineLayer;
using MassSpectrometry;
using Microsoft.Win32;
using Newtonsoft.Json;
using Proteomics;
using Readers;
using System.Collections.Generic;
using System.Data;
using System.IO;
using System.Linq;
using System.Windows;
using TaskLayer;
using TorchSharp;
using UsefulProteomicsDatabases;

namespace MultiModSearchViz
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {

        public List<DataTable> RunResults { get; set; }
        public MainWindow()
        {
            InitializeComponent();
        }

        private void LoadPsm(object sender, RoutedEventArgs e)
        {
            OpenFileDialog psmFile = new OpenFileDialog();
            psmFile.ShowDialog();
            psmFile.FileName = @"Asdad";

            //var psmList = MultiModSearch.ReadFilteredPsmTSVShort(psmFile.FileName);

        }

        private void LoadDataFile(object sender, RoutedEventArgs e)
        {
            OpenFileDialog dataFilePath = new OpenFileDialog();
            dataFilePath.ShowDialog();
            dataFilePath.FileName = @"asdada";

            var dataFile = MsDataFileReader.GetDataFile(dataFilePath.FileName).GetAllScansList()
                .Where(x => x.MsnOrder == 2).ToArray();
        }

        private void RunButton(object sender, RoutedEventArgs e)
        {

            List<FilteredPsmTSV> psmList = MultipleSearchEngine.ReadFilteredPsmTSVShort($@"{psmPath.Text}");

            MsDataFile dataFile = MsDataFileReader.GetDataFile($@"{mzMLPath.Text}");

            GlobalVariables.SetUpGlobalVariables();

            int maxNumOfMods = int.Parse(this.maxNumOfMods.Text);

            List<Modification> commonBiologycalMods = GlobalVariables.AllModsKnown.OfType<Modification>()
                .Where(mod => mod.ModificationType.Contains("Common Biological") ||
                              mod.ModificationType.Contains("Common Artifact") ||
                              mod.ModificationType.Contains("Common Variable")).ToList();



            List<Modification> modsFromUnimod =
                Loaders.LoadUnimod(@"../../../../EngineLayer/unimod.xml").ToList();
            //Loaders.LoadUnimod(
            //        @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
            //    .ToList();

            List<Modification> fixedMods = new List<Modification>();
            fixedMods.Add(modsFromUnimod.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")));

            new MultipleSearchEngine(psmList, commonBiologycalMods, maxNumOfMods,
                fixedMods, dataFile, true).Run(pathToSave.Text, fixedMods, dataFile);

            //var (tables, objects) = MultipleSearchEngine.Run(engine, psmList,
            //    fixedMods, dataFile, maxNumOfMods, pathToSave.Text);

            //dataGrid.ItemsSource = tables.Find(x => x.TableName.Equals(proteinGroups.SelectedItems)).Rows;
            //proteinGroups.ItemsSource = run.Select(x => x.TableName).ToList();
            //dataGrid.DataContext = run;
            //RunResults = tables;
        }

        private void VizResults(object sender, RoutedEventArgs e)
        {
            string jsonString = File.ReadAllText(pathToResults.Text);
            
            List<MultipleSearchResults> resultsFile = System.Text.Json.JsonSerializer
                .Deserialize<List<MultipleSearchResults>>(jsonString);

            resultsFile.RemoveAll(x => x == null);
            
            var groupedPeptides = resultsFile.GroupBy(x => x.BaseSequence).ToList();
            
            RunResults = MultipleSearchResults.GetDataTables(groupedPeptides);
            
            proteinGroups.ItemsSource = groupedPeptides.Select(x => x.Key).ToList();
        }

        private void UpdatePeptideGroup(object sender, RoutedEventArgs e)
        {

            var peptideGroup = proteinGroups.SelectedItem.ToString();

            var dataView = RunResults.Find(x => x.TableName.Equals(peptideGroup)).AsDataView();
            dataGrid.ItemsSource = dataView;
        }
    }
}
