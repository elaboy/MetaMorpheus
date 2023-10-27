using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CsvHelper;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using NUnit.Framework;
using Org.BouncyCastle.Asn1.Cmp;
using TorchSharp;

namespace Test
{
    public sealed class TestTensorize
    {
        [Test]
        public void TestTensorLoadingFromCSVChronologerDB()
        {

            var configuration = new CsvConfiguration(CultureInfo.InvariantCulture)
            {
                HasHeaderRecord = true,
                Delimiter = "\t"
            };


            var reader =
                new StreamReader(@"F:\Research\Data\Chronologer\Chronologer_DB_220308\Chronologer_DB_220308.csv");

            var csv = new CsvReader(reader, configuration);

            var records = csv.GetRecords<ChronologerDB>().ToList();

            reader.Close();
            csv.Dispose();

        }

        [Test]
        public void TestTensorizeDB()
        {
            var db = ChronologerDB.LoadChronologerDB();

            List<torch.Tensor> tensorList = new List<torch.Tensor>();

            for (int i = 0; i < 2; i++)
            {
                tensorList.Add(torch.zeros(2, db[i].PeptideModSeq.Length));
            }

            Console.Write(tensorList[0].ToString(TensorStringStyle.Numpy));

        }
    }

    public class ChronologerDB
    {
        [Index(0)]
        public string Source { get; set; }
        [Index(1)]
        public string PeptideModSeq { get; set; }
        [Index(2)]
        public double Prosit_RT { get; set; }
        [Index(3)]
        public double HI { get; set; }

        public static List<ChronologerDB> LoadChronologerDB()
        {
            var configuration = new CsvConfiguration(CultureInfo.InvariantCulture)
            {
                HasHeaderRecord = true,
                Delimiter = "\t"
            };

            var reader =
                new StreamReader(@"F:\Research\Data\Chronologer\Chronologer_DB_220308\Chronologer_DB_220308.csv");

            var csv = new CsvReader(reader, configuration);

            var records = csv.GetRecords<ChronologerDB>().ToList();

            reader.Close();
            csv.Dispose();

            return records;
        }

        public static List<torch.Tensor> Tensorize(List<ChronologerDB> db)
        {

            List<torch.Tensor> tensorList = new List<torch.Tensor>();

            for (int i = 0; i < 2; i++)
            {
                tensorList.Add(torch.zeros(2, db[i].PeptideModSeq.Length));
            }

            return tensorList;
        }
    }
    
}


