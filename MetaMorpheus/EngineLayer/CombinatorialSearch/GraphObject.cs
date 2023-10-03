using System.Collections.Generic;
using Proteomics;

namespace EngineLayer.CombinatorialSearch
{
    /// <summary>
    /// Network with all the Nodes (Residues)
    /// </summary>
    public class GraphObject
    {
        public List<Node> Nodes { get; set; }
        public string BaseSequence { get; set; }
        public double Coverage { get; set; }
        public List<List<Node>> ConnectedNodes { get; set; } //new list for each fragment
        public GraphObject(PeptideSpectralMatch psm)
        {
            Nodes = new List<Node>();
            foreach (var residue in psm.BaseSequence)
            {
                Nodes.Add(new Node(residue));
            }
        }
    }

    /// <summary>
    /// Each Node represents each residue
    /// </summary>
    public class Node
    {
        public char Residue { get; set; }
        public Modification Modification { get; set; }
        public double? MZ { get; set; }
        public bool Matched { get; set; }
        public bool DatabaseEvidence { get; set; }
        public Node(char residue)
        {
            Matched = false;
            DatabaseEvidence = false;
        }
    }
}
