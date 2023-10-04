using System.Collections.Generic;
using System.Collections.Immutable;
using System.Dynamic;
using System.Linq;
using System.Reflection;
using MassSpectrometry;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.CombinatorialSearch
{
    /// <summary>
    /// Network with all the Nodes (Residues)
    /// </summary>
    public class GraphObject
    {
        public ImmutableList<Node> Nodes { get; set; }
        public PeptideSpectralMatch PSM { get; set; }
        public PeptideWithSetModifications Peptide { get; set; }
        public List<Product> Products { get; set; }
        public List<Product> NterminalProducts { get; set; }
        public List<Product> CterminalProducts { get; set; }
        public string BaseSequence { get; set; }
        public Dictionary<int, Modification> Modifications { get; set; }
        public int CterminusIndex { get; set; }
        public int NterminusIndex { get; set; }
        public double Coverage { get; set; }
        public List<List<Node>> ConnectedNodes { get; set; } //new list for each fragment

        public GraphObject(PeptideSpectralMatch psm)
        {
            PSM = psm;
            var nodes = new List<Node>();

            for (int i = 0; i < psm.BaseSequence.Length; i++)
            {
                nodes.Add(new Node(psm.BaseSequence[i], i+1,
                    psm.BaseSequence.Length-(i)));
            }

            Nodes = nodes.ToImmutableList();
            Peptide = psm.BestMatchingPeptides.First().Peptide;
            Products = new List<Product>();
            Modifications = Peptide.AllModsOneIsNterminus;

            //Sort Products into separate lists
            NterminalProducts = new List<Product>();
            CterminalProducts = new List<Product>();
            FragmentAndSortNodeProducts();
            SetMods();

            CterminusIndex = 0;
            NterminusIndex = Nodes.Count;
        }

        private void FragmentAndSortNodeProducts()
        {
            Peptide.Fragment(DissociationType.HCD, FragmentationTerminus.N, NterminalProducts);
            Peptide.Fragment(DissociationType.HCD, FragmentationTerminus.C, CterminalProducts);

            for (int i = 0; i < NterminalProducts.Count; i++)
            {
                Nodes[i].NterminusTheoreticalProduct = NterminalProducts[i];
            }

            if (Nodes.Last().Residue.Equals('L') || Nodes.Last().Residue.Equals('L'))
            {
                for (int i = CterminalProducts.Count - 1; i > 1; i--)
                {
                    Nodes[i].CterminusTheoreticalProduct = CterminalProducts[i];
                }
            }
            else
            {
                for (int i = CterminalProducts.Count; i > 1; i--)
                {
                    Nodes[i].CterminusTheoreticalProduct = CterminalProducts[i];
                }
            }
        }

        private void SetMods()
        {
            foreach (var mod in Modifications)
            {
                Nodes[mod.Key].Modification = mod.Value;
            }
        }
        public void MoveInwards()
        {
            CterminusIndex = CterminusIndex - 1;
            NterminusIndex = NterminusIndex + 1;
        }

        public bool CheckEdgeNodes()
        {
            if (Nodes[CterminusIndex].Matched && Nodes[NterminusIndex].Matched)
                return true;

            return false;
        }
    }

    /// <summary>
    /// Each Node represents each residue
    /// </summary>
    public class Node
    {
        public char Residue { get; set; }
        public int NterminusIndex { get; set; }
        public int CterminusIndex { get; set; }
        public Product CterminusTheoreticalProduct { get; set; }
        public Product NterminusTheoreticalProduct { get; set; }
        public Product CterminusExperimentalProduct { get; set; }
        public Product? NterminusExperimentalProduct { get; set; }
        public Modification Modification { get; set; }
        public double? EvidencedMZ { get; set; }
        public bool Matched { get; set; }
        public bool DatabaseEvidence { get; set; }
        public Node(char residue, int nterminus, int cterminus)
        {
            Residue = residue;
            Matched = false;
            DatabaseEvidence = false;
            NterminusIndex = nterminus;
            CterminusIndex = cterminus;
        }
    }
}
