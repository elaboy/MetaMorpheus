using EngineLayer;
using MassSpectrometry;
using Microsoft.Win32;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Threading.Tasks;
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
                .Where(mod => mod.ModificationType.Contains("Common Biological")).ToList();

            List<Modification> modsFromUnimod =
                Loaders.LoadUnimod(
                        @"C:\Users\Edwin\Documents\GitHub\MetaMorpheus\MetaMorpheus\EngineLayer\Data\unimod.xml")
                    .ToList();

            List<Modification> fixedMods = new List<Modification>();
            fixedMods.Add(modsFromUnimod.Find(x => x.IdWithMotif.Equals("Carbamidomethyl on C")));

            var engine = new MultipleSearchEngine(commonBiologycalMods, maxNumOfMods, true);

            var run = MultipleSearchEngine.Run(engine, psmList, fixedMods, dataFile, maxNumOfMods);

            proteinGroups.ItemsSource = run.Select(x => x.TableName).ToList();
            dataGrid.ItemsSource = run.Select(x => x.Rows[0]);
        }
    }
}
