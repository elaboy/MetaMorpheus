using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using Microsoft.ML.Trainers.FastTree;
using Proteomics.Fragmentation;
using Easy.Common.Extensions;

namespace EngineLayer.CombinatorialSearch
{
    public static class CSExtension
    {

        public static void BuildCombinationWithAddedMass(List<List<Modification>> listOfModCombination,
            out List<KeyValuePair<double, Modification[]>> combinationsWithAddedMass)
        {
            combinationsWithAddedMass = listOfModCombination
                .Select(x => new KeyValuePair<double, Modification[]>(
                key: x.Select(x => x.MonoisotopicMass.Value).Sum(),
                value: x.Select(x => x).ToArray()))
                .OrderBy(x => x.Key).ToList();
        }

        /// <summary>
        /// Recursive method that builds the combination of mods. 
        /// </summary>
        /// <param name="sortedListOfModsToAdd"></param>
        /// <param name="listOfModCombinations"></param>
        /// <param name="maxNumberOfModsInGroup"></param>
        /// <param name="allModsFromOneToN"></param>
        public static void CombinationBuilder(List<Modification> sortedListOfModsToAdd,
            ref List<List<Modification>> listOfModCombinations, int maxNumberOfModsInGroup, bool allModsFromOneToN)
        {
            List<List<Modification>> newModLists = new();
            if (maxNumberOfModsInGroup == 0)
            {
                foreach (var modificationToAdd in sortedListOfModsToAdd)
                {
                    newModLists.Add(new List<Modification>() { modificationToAdd });
                }
                newModLists = newModLists.DistinctBy(n => ModListNameString(n)).ToList();
                listOfModCombinations = newModLists;
            }
            else
            {
                CombinationBuilder(sortedListOfModsToAdd, ref listOfModCombinations, maxNumberOfModsInGroup - 1, allModsFromOneToN);
                newModLists.Clear();
                foreach (var modList in listOfModCombinations.Where(c => c.Count == (maxNumberOfModsInGroup - 1)))
                {
                    foreach (var modificationToAdd in sortedListOfModsToAdd)
                    {
                        List<Modification> newModList = modList.ToList();
                        newModList.Add(modificationToAdd);
                        newModLists.Add(newModList.OrderBy(n => n.IdWithMotif).ToList());
                    }
                }
                newModLists = newModLists.DistinctBy(n => ModListNameString(n)).ToList();
                listOfModCombinations.AddRange(newModLists);
                if (!allModsFromOneToN)
                {
                    listOfModCombinations = listOfModCombinations.Where(c => c.Count == maxNumberOfModsInGroup).ToList();
                }
            }

        }
        /// <summary>
        /// Returns the Modifications IdWithMotif joined. For CombinationBuilder use.
        /// </summary>
        /// <param name="list"></param>
        /// <returns></returns>
        private static string ModListNameString(List<Modification> list)
        {
            return String.Join("", list.Select(n => n.IdWithMotif));
        }

        //For GraphObject Use only
        public static void ResetNodeMatchBoolToFalse(GraphObject network)
        {
            foreach (var node in network.Nodes)
            {
                node.CterminusMatched = false;
                node.NterminusMatched = false;
            }
        }

        public static void SetMods(GraphObject network, Dictionary<int, Modification> mods)
        {
            foreach (var mod in mods)
            {
                // network.Nodes[mod.Key].Modification = mod.Value;
            }
        }
        public static void FragmentAndSortNodeProducts(GraphObject network)
        {
            network.Peptide.Fragment(DissociationType.HCD, FragmentationTerminus.N, network.NterminalProducts);
            network.Peptide.Fragment(DissociationType.HCD, FragmentationTerminus.C, network.CterminalProducts);
            int breasss = 0;
            for (int i = 0; i < network.NterminalProducts.Count; i++)
            {
                if(i ==0)
                    network.Nodes[i + 1].NterminusTheoreticalProduct = network.NterminalProducts[i];
                else
                    network.Nodes[i].NterminusTheoreticalProduct = network.NterminalProducts[i];
            }

            network.CterminalProducts.Reverse();
            for (int i = 0; i < network.CterminalProducts.Count; i++)
            {
                if (i == 0)
                    network.Nodes[i + 1].CterminusTheoreticalProduct = network.CterminalProducts[i];
                else
                    network.Nodes[i].CterminusTheoreticalProduct = network.CterminalProducts[i];
            }
        }

        public static void FixMissedCleavageOnTerminal(GraphObject network)
        {
            if (network.Nodes.Last().CterminusTheoreticalProduct.Annotation.Equals("a0"))
            {
                foreach (var node in network.Nodes)
                {
                    node.CterminusIndex = node.CterminusIndex - 1;
                }
                MatchSpectraAndScore(network);
            }
        }
        public static void MatchSpectraAndScore(GraphObject network)
        {
            //Nterminal
            var NterminalMatch = MetaMorpheusEngine.MatchFragmentIons(
                new Ms2ScanWithSpecificMass(network.PSM.MsDataScan,
                    network.PSM.ScanPrecursorMonoisotopicPeakMz,
                    network.PSM.ScanPrecursorCharge, network.PSM.FullFilePath, new CommonParameters()),
                network.NterminalProducts,
                new CommonParameters(), true);

            //Cterminal
            network.CterminalProducts.Reverse(); //reverts the products order to begin with 1
            var CterminalMatch = MetaMorpheusEngine.MatchFragmentIons(
                new Ms2ScanWithSpecificMass(network.PSM.MsDataScan,
                    network.PSM.ScanPrecursorMonoisotopicPeakMz,
                    network.PSM.ScanPrecursorCharge, network.PSM.FullFilePath, new CommonParameters()),
                network.CterminalProducts,
                new CommonParameters(), true);
            network.CterminalProducts.Reverse(); //back to reverse

            //Update matched status and sort MatchedFragmentIon to each Node
            foreach (var node in network.Nodes)
            {
                node.NterminusFragmentIon = new List<MatchedFragmentIon>();
                foreach (var fragment in NterminalMatch)
                {
                    if (node.NterminusIndex.Equals(fragment.NeutralTheoreticalProduct.FragmentNumber))
                    {
                        node.NterminusFragmentIon.Add(fragment);
                        node.NterminusMatched = true;
                        break;
                    }
                }
            }

            foreach (var node in network.Nodes)
            {
                node.CterminusFragmentIon = new List<MatchedFragmentIon>();
                foreach (var fragment in CterminalMatch)
                {
                    if (node.CterminusIndex.Equals(fragment.NeutralTheoreticalProduct.FragmentNumber))
                    {
                        node.CterminusFragmentIon.Add(fragment);
                        node.CterminusMatched = true;
                        break;
                    }
                }
            }

            //Update Coverage todo
            network.NetworkCoverage = 0; //reset to 0

            foreach (var node in network.Nodes)
            {
                if (node.CterminusMatched)
                    network.NetworkCoverage = network.NetworkCoverage + 1;
                if (node.NterminusMatched)
                    network.NetworkCoverage = network.NetworkCoverage + 1;
            }

            network.NetworkCoverage = network.NetworkCoverage / (network.NterminalProducts.Count + network.CterminalProducts.Count);
        }
        private static IEnumerable<Dictionary<int, Modification>> GetModDictionary(GraphObject network, List<Modification> modifications)
        {
            var listOfModsDict = new List<Dictionary<int, Modification>>();

            foreach (var mod in network.ModsToTry)
            {
                var moddedPeptide = network.Peptide.Protein
                    .Digest(network.Peptide.DigestionParams,
                        new List<Modification>(), mod);

                var peptides = moddedPeptide.Select(x =>
                    x.BaseSequence.Equals(network.BaseSequence) ? x : null);

                var trash = new List<byte>();


                foreach (var pep in peptides)
                {
                    if (pep is null) continue;

                    listOfModsDict.Add(pep.AllModsOneIsNterminus);
                }
                // peptides.ForEach(x => listOfModsDict.Add(x.AllModsOneIsNterminus));
            }

            return listOfModsDict.Where(x => x.IsNotNullOrEmpty());
        }

    }
}
