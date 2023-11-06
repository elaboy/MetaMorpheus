using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteomics;
using Proteomics.Fragmentation;

namespace EngineLayer.CombinatorialSearch
{
    public class NodeNet
    {
        public PeptideSpectralMatch Psm { get; set; }
        public List<Tree> Trees { get; set; }

        public NodeNet()
        {

        }
    }

    public class Tree
    {
        public List<Leaf> Leaves { get; set; }
        public List<Modification> TreeModifications { get; set; }
        public Tree()
        {

        }

    }

    public class Leaf
    {
        public int LeafLevel { get; set; }
        public List<Product> Products { get; set; }

        public Leaf()
        {

        }
    }
}
