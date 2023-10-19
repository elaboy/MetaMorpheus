using EngineLayer;
using EngineLayer.CombinatorialSearch;
using NUnit.Framework;
using System.Collections.Generic;
using TaskLayer;
using TaskLayer.CombinatorialSearchTask;

namespace Test
{
    public sealed class TestNodeNetwork
    {
        public const string TestingFile = @"F:\Research\Data\Search_Algorithm\OneMs1Files\test.mzML";
        [Test]
        public void TestNodeNetworkConstructor()
        {
            var taskList = new List<(string, MetaMorpheusTask)>();
            var CSTask = new CSTask(new CommonParameters());

            taskList.Add(("CS-Task", CSTask));

            List<DbForTask> dbForTask = new List<DbForTask>();
            dbForTask.Add(new DbForTask(
                @"F:\Research\Data\Search_Algorithm\08-30-22_bottomup\uniprotkb_taxonomy_id_9606_AND_reviewed_2023_09_04.fasta",
                false));

            var runner = new EverythingRunnerEngine(taskList,
                new List<string>() { TestingFile},
                dbForTask, @"F:\TestingCSTask");

            runner.Run();
        }
    }
}
