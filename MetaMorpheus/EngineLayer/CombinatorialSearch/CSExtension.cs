using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

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
    }
}
