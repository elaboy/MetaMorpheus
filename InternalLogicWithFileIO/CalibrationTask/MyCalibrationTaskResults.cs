﻿using System.Text;

namespace InternalLogicTaskLayer
{
	class MyCalibrationTaskResults : MyTaskResults
	{
		public MyCalibrationTaskResults(MyTaskEngine s) : base(s)
		{
		}

		protected override string GetStringForOutput()
		{
			var sb = new StringBuilder();
			return sb.ToString();
		}
	}
}