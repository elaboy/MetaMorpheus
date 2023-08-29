using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;

namespace TaskLayer.CombinatorialSearchTask
{
    public class CSTask : MetaMorpheusTask
    {
        public CommonParameters CommonParameters { get; set; }
        public CSTask() : base(MyTask.CombinatorialSearch)
        {

        }

        protected override MyTaskResults RunSpecific(string OutputFolder,
            List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId,
            FileSpecificParameters[] fileSettingsList)
        {
            // setup 

            // search 

            // CS stuff 


            // output results


            throw new NotImplementedException();
        }
    }
}
