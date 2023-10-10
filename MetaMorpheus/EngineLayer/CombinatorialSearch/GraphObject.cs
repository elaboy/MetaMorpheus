﻿using MassSpectrometry;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Diagnostics;
using System.Linq;

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
        public List<Product> NterminalProducts { get; set; }
        public List<Product> CterminalProducts { get; set; }
        public string BaseSequence { get; set; }
        public Dictionary<int, Modification> Modifications { get; set; }
        public int NetworkCterminusIndex { get; set; }
        public int NetworkNterminusIndex { get; set; }
        public double NetworkCoverage { get; set; }
        public List<List<Node>> ConnectedNodes { get; set; } //new list for each fragment

        public GraphObject(PeptideSpectralMatch psm) //todo pass common parameters to feed the experiment specifics such as DissociationType
        {
            PSM = psm;
            var nodes = new List<Node>();

            for (int i = 0; i < psm.BaseSequence.Length; i++)
            {
                nodes.Add(new Node(psm.BaseSequence[i], i + 1,
                    psm.BaseSequence.Length - (i)));
            }

            Nodes = nodes.ToImmutableList();
            // psm.BestMatchingPeptides.First().Peptide.AllModsOneIsNterminus.Clear(); //no mods
            Peptide = psm.BestMatchingPeptides.First().Peptide;
            Modifications = Peptide.AllModsOneIsNterminus;

            //Sort Products into separate lists
            NterminalProducts = new List<Product>();
            CterminalProducts = new List<Product>();
            FragmentAndSortNodeProducts();
            SetMods();

            NetworkNterminusIndex = 0;
            NetworkCterminusIndex = Nodes.Count;
            NetworkCoverage = 0;

            MatchSpectraAndScore();
            FixMissedCleavageOnTerminal();
            CheckEdgeNodes();
        }

        private void SetModificationDatabase()
        {

        }

        private void MatchSpectraAndScore()
        {
            //Nterminal
            var NterminalMatch = MetaMorpheusEngine.MatchFragmentIons(
                new Ms2ScanWithSpecificMass(PSM.MsDataScan,
                    PSM.ScanPrecursorMonoisotopicPeakMz,
                    PSM.ScanPrecursorCharge, PSM.FullFilePath, new CommonParameters()),
                NterminalProducts,
                new CommonParameters(), true);

            //Cterminal
            CterminalProducts.Reverse(); //reverts the products order to begin with 1
            var CterminalMatch = MetaMorpheusEngine.MatchFragmentIons(
                new Ms2ScanWithSpecificMass(PSM.MsDataScan,
                    PSM.ScanPrecursorMonoisotopicPeakMz,
                    PSM.ScanPrecursorCharge, PSM.FullFilePath, new CommonParameters()),
                CterminalProducts,
                new CommonParameters(), true);

            //Update matched status and sort MatchedFragmentIon to each Node
            foreach (var node in Nodes)
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

            foreach (var node in Nodes)
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
            NetworkCoverage = 0; //reset to 0

            foreach (var node in Nodes)
            {
                if (node.CterminusMatched)
                    NetworkCoverage = NetworkCoverage + 1;
                if(node.NterminusMatched)
                    NetworkCoverage = NetworkCoverage + 1;
            }

            NetworkCoverage = NetworkCoverage / (NterminalProducts.Count + CterminalProducts.Count);
        }

        private void FixMissedCleavageOnTerminal()
        {
            if (Nodes.Last().CterminusTheoreticalProduct.Annotation.Equals("a0"))
            {
                foreach (var node in Nodes)
                {
                    node.CterminusIndex = node.CterminusIndex - 1;
                }
                MatchSpectraAndScore();
            }

        }

        private void FragmentAndSortNodeProducts()
        {
            Peptide.Fragment(DissociationType.HCD, FragmentationTerminus.N, NterminalProducts);
            Peptide.Fragment(DissociationType.HCD, FragmentationTerminus.C, CterminalProducts);

            for (int i = 0; i < NterminalProducts.Count; i++)
            {
                Nodes[i].NterminusTheoreticalProduct = NterminalProducts[i];
            }

            CterminalProducts.Reverse();
            for (int i = 0; i < CterminalProducts.Count; i++)
            {
                Nodes[i].CterminusTheoreticalProduct = CterminalProducts[i];
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
            NetworkCterminusIndex = NetworkCterminusIndex - 1;
            NetworkNterminusIndex = NetworkNterminusIndex + 1;
        }

        public bool CheckEdgeNodes() //todo fix for annotation check
        {
            if ((Nodes[NetworkNterminusIndex].NterminusMatched ||
                 Nodes[NetworkNterminusIndex].NterminusFragmentIon.Equals("a0")) &&
                (Nodes[NetworkNterminusIndex].CterminusMatched ||
                 Nodes[NetworkNterminusIndex].CterminusFragmentIon.Equals("a0")) &&
                (Nodes[NetworkCterminusIndex].NterminusMatched ||
                 Nodes[NetworkCterminusIndex].NterminusFragmentIon.Equals("a0")) &&
                (Nodes[NetworkCterminusIndex].CterminusMatched ||
                 Nodes[NetworkCterminusIndex].CterminusFragmentIon.Equals("a0")))

            {
                Debug.WriteLine("Edges matched, moving inward");
                return true;

            }

            if (Nodes[NetworkCterminusIndex - 1].CterminusTheoreticalProduct.NeutralMass >
                PSM.MsDataScan.MassSpectrum.LastX)
            {
                Debug.WriteLine("mz over limit, moving inward");
                return true;
            }

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
        public List<MatchedFragmentIon> NterminusFragmentIon { get; set; }
        public List<MatchedFragmentIon> CterminusFragmentIon { get; set; }
        public Modification Modification { get; set; }
        public double? EvidencedMZ { get; set; }
        public bool NterminusMatched { get; set; }
        public bool CterminusMatched { get; set; }
        public bool DatabaseEvidence { get; set; }
        public Node(char residue, int nterminus, int cterminus)
        {
            Residue = residue;
            NterminusMatched = false;
            CterminusMatched = false;
            DatabaseEvidence = false;
            NterminusIndex = nterminus;
            CterminusIndex = cterminus;
        }
    }
}
