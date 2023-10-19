using System;
using Proteomics;
using System.Collections.Generic;
using System.Dynamic;
using System.Linq;
using Easy.Common.Extensions;
using FlashLFQ;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.CombinatorialSearch
{
    public class NodeNetwork
    {
        public List<PeptideWithSetModifications> PeptidesWithSetMods { get; set; }
        public List<List<Modification>> NetworkModifications{ get; private set; }
        public List<Dictionary<int, Modification>> ModsForNodes { get; private set; }
        public static IEnumerable<IGrouping<int, Dictionary<int, Modification>>> GroupedModsByFirstMod { get; private set; } 
        public IEnumerable<IGrouping<int, NodeStructure>> NodesByResidue { get; private set; }
        public List<NodeStructure> ResidueNodes { get; private set; }
        public static PeptideSpectralMatch PSM { get; private set; }
        public NodeNetwork(PeptideSpectralMatch psm, List<List<Modification>> modsThatFitDelta)
        {
            PSM = psm;
            NetworkModifications = modsThatFitDelta;
            PeptidesWithSetMods = new List<PeptideWithSetModifications>();
            ModsForNodes = new List<Dictionary<int, Modification>>();
            ResidueNodes = new List<NodeStructure>();
            SetUpNet();
        }

        private void SetUpNet()
        {
            foreach (var mod in NetworkModifications)
            {
                var peptideAsProtein = new Protein(PSM.BaseSequence, PSM.ProteinAccession, PSM.Organism)
                    .Digest(new DigestionParams("top-down", 0), new List<Modification>(),
                        mod);

                foreach (var variablePeptide in peptideAsProtein)
                {
                    PeptidesWithSetMods.Add(variablePeptide);
                }
            }

            PeptidesWithSetMods = PeptidesWithSetMods.DistinctBy(x => x.FullSequence).ToList();
            PeptidesWithSetMods.RemoveAt(0); //remove first empty dict
            ExtractModificationDictionary();

            //Make Parent Node
            foreach (var group in GroupedModsByFirstMod)
            {
                for (int i = 0; i < PSM.BaseSequence.Length; i++)
                {
                    foreach (var variable in group)
                    {
                        ResidueNodes.Add(new NodeStructure(PSM.BaseSequence[i], i, variable));
                    }
                }
            }

            NodesByResidue = ResidueNodes.GroupBy(x => x.ResiduePosition);

            //make children nodes
            foreach (var node in NodesByResidue)
            {
                foreach (var variation in node)
                {
                    int level = 1;
                    RecursivelyAddChildren(variation, level);
                }
            }
        }

        private void RecursivelyAddChildren(NodeStructure variation, int level)
        {
            int i = level;
            while (i < PSM.BaseSequence.Length + 1)
            {
                if (variation.SubSequence.Length > 1)
                {
                    if (variation.ModificationsForThisTree.ContainsKey(i))
                    {
                        //take parent mods and add this position's mod
                        Dictionary<int, Modification> modsForThisLevel = variation.Parent.ModificationsForThisNode;
                        modsForThisLevel.Add(i, variation.Parent.ModificationsForThisTree[i]);

                        //p^T-1 = number of nodes to be generated
                        var numberOfModdableResidues = 0;
                        foreach (var residue in variation.SubSequence)
                        {
                            foreach (var mod in modsForThisLevel)
                            {
                                if (char.Parse(mod.Value.Target.ToString()) == residue)
                                {
                                    numberOfModdableResidues = numberOfModdableResidues + 1;
                                }
                            }
                        }

                        if (numberOfModdableResidues <= 1)
                        {
                            new NodeStructure(variation, variation.ModificationsForThisTree,
                                modsForThisLevel, i);
                            i = i + 1;
                            if (i >= PSM.BaseSequence.Length+1) break;

                            RecursivelyAddChildren(variation.Children[0], i);
                        }

                        //p^T-1
                        int numberOfNodesToBeGenerated =
                            (numberOfModdableResidues ^
                             (modsForThisLevel.GroupBy(x => x.Value.Target).Distinct().Count())) - 1;

                        //generate nodes for this level
                        for (int k = 0; k < numberOfNodesToBeGenerated; k++)
                        {
                            new NodeStructure(variation, variation.ModificationsForThisTree,
                                modsForThisLevel, i);
                        }
                    }
                    else
                    {
                        new NodeStructure(variation, variation.ModificationsForThisTree,
                            variation.ModificationsForThisNode, i);
                    }
                }
                else
                {
                    //continue;
                    new NodeStructure(variation, variation.ModificationsForThisTree,
                        new Dictionary<int, Modification>(), i);
                }
                i = i + 1;
                if (i >= PSM.BaseSequence.Length+1) break;

                RecursivelyAddChildren(variation.Children[0], i);
            }
        }



        private void ExtractModificationDictionary()
        {
            foreach (var peptide in PeptidesWithSetMods)
                ModsForNodes.Add(peptide.AllModsOneIsNterminus);

            GroupedModsByFirstMod = ModsForNodes.GroupBy(x => x.Keys.First()).OrderBy(x => x.Key);
        }
    }

    public class NodeStructure
    {
        public NodeStructure? Parent { get; private set; }
        public Dictionary<int, Modification> ModificationsForThisTree { get; private set; }
        public Dictionary<int, Modification> ModificationsForThisNode { get; private set; }
        public string SubSequence { get; private set; }
        public char? Residue { get; private set; }
        public int ResiduePosition { get; set; } //0 means N-terminus
        public List<NodeStructure> Children { get; private set; }
        public int? Level { get; private set; }

        //parent
        public NodeStructure(char? residue, int residuePosition, Dictionary<int, Modification> modsForThisLevel)
        {
            Residue = residue;
            ResiduePosition = residuePosition;
            ModificationsForThisTree = modsForThisLevel;
            Level = 0;
            SubSequence = Residue.ToString();
        }
        //children
        public NodeStructure(NodeStructure parentNode, Dictionary<int, Modification> modDictionaryForNewParent,
            Dictionary<int, Modification> modsForThisLevel, int treeLevel)
        {
            parentNode.Children = new List<NodeStructure>();
            Parent = parentNode;
            Level = treeLevel;
            ModificationsForThisNode = modsForThisLevel;
            ModificationsForThisTree = modDictionaryForNewParent;
            var subsequence = NodeNetwork.PSM.BaseSequence.Substring(0, treeLevel);
            SubSequence = subsequence;
            parentNode.Children.Add(this);
            //var targets = new List<ModificationMotif>();

            //foreach (var mod in modDictionaryForNewParent)
            //{
            //    if (!targets.Contains(mod.Value.Target))
            //        targets.Add(mod.Value.Target);
            //}

            //for (int i = 0; i < NodeNetwork.PSM.BaseSequence.Length; i++)
            //{



            //}


            //int nodesQuant = 0;

            //foreach (var target in targets)
            //{
            //    foreach (var residue in NodeNetwork.PSM.BaseSequence)
            //        if (residue == char.Parse(target.ToString()))
            //            nodesQuant = nodesQuant + 1;
            //}

            //MakeChildren(parentNode);
        }

        public NodeStructure() { }

        //private void MakeChildren(NodeStructure parentNode, int numberOfNodes)
        //{
        //    parentNode.Children = new NodeStructure[numberOfNodes];
        //}
    }

}



//public class Map
//{
//    public static MsDataScan Scan { get; set; }
//    private List<SwappablePeptideWithSetModifications> peptides { get; set; }
//    private List<ParentNode> parentNodes { get; set; }
//    private Dictionary<int, (bool, bool)> fragmentMatch { get; set; }
//    public Map(PeptideSpectralMatch psm, List<List<Modification>> modsThatFitDelta, List<Modification> variableMods = null)
//    {
//        if(variableMods is null)
//            variableMods = new List<Modification>();
//        Scan = psm.MsDataScan;
//        peptides = new List<SwappablePeptideWithSetModifications>();
//        var listOfProducts = new List<Product>();

//        //Dictionary with boolean for each fragment terminal (bool, bool) (Nterminus, Cterminus)
//        fragmentMatch = new Dictionary<int, (bool, bool)>();
//        for (int i = 0; i < psm.BaseSequence.Length; i++)
//        {
//            fragmentMatch.Add(i, (false,false));
//        }

//        foreach (var mod in modsThatFitDelta)
//        {
//            var peptideAsProtein = new Protein(psm.BaseSequence, psm.ProteinAccession, psm.Organism).Digest(
//                new DigestionParams("top-down", 0), mod, 
//                variableMods);

//            foreach (var variablePeptide in peptideAsProtein)
//            {
//                peptides.Add(new SwappablePeptideWithSetModifications(variablePeptide, variablePeptide.AllModsOneIsNterminus));
//            }
//        }

//        var modifiedIndecies = peptides.SelectMany(p => p.AllModsOneIsNterminus.Keys).Distinct().OrderBy(p => p).ToArray();
//        var indexedNodes = new List<Dictionary<int, Modification>>();
//        foreach (var index in modifiedIndecies)
//        {
//            foreach (var peptide in peptides)
//            {
//                if (peptide.AllModsOneIsNterminus.ContainsKey(index))
//                {
//                    indexedNodes.Add(peptide.AllModsOneIsNterminus);
//                }
//            }
//        }



//        foreach (var groupByIndex in peptides.GroupBy(p => p.AllModsOneIsNterminus.Keys.Min()))
//        {

//        }



//        // var groupedByMods = peptides.GroupBy(x => x.AllModsOneIsNterminus.Keys.Min());
//        var orderedIndexedNodes = indexedNodes.GroupBy(x => x.Keys.Min());
//        parentNodes = new List<ParentNode>();

//        //foreach (var group in groupedByMods)
//        //{
//        //    parentNodes.Add(new ParentNode(group));
//        //}
//    }

//    private void DummyRecursion(List<SwappablePeptideWithSetModifications> list, int[] modifiedIndexes,
//        int currentIndex)
//    {
//        if (currentIndex > modifiedIndexes.Length || !list.Any())
//            return;


//        int currentModIndex = modifiedIndexes[currentIndex];
//        foreach (var modGroup in list.Where(p => p.AllModsOneIsNterminus.ContainsKey(currentModIndex))
//                     .GroupBy(p => p.AllModsOneIsNterminus[currentModIndex].IdWithMotif))
//        {
//            // construct node per group

//            DummyRecursion(modGroup.ToList(), modifiedIndexes, currentIndex + 1);
//        }
//    }
//}