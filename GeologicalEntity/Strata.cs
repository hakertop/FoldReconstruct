using OSGeo.OGR;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeologicalEntity
{
    public class Strata
    {
        /// <summary>
        /// All strata
        /// </summary>
        public List<Stratum> pAllStrata { get; set; }

        public Strata(Dictionary<int, Stratum> _pAllStratums)
        {
            pAllStrata = new List<Stratum>();
            foreach (var vs in _pAllStratums.Values)
                pAllStrata.Add(vs);

            pAllStrata.TrimExcess();


        }

        /// <summary>
        /// Topology relation database
        /// </summary>
        /// <param name="pAllStrata"></param>
        /// <returns></returns>
        public Dictionary<string,string> CreateTopoDatabase()
        {
            Dictionary<string, string> pStrataTopoDatabase = new Dictionary<string, string>();

            foreach (var s1 in pAllStrata)
            {
                foreach(var s2 in pAllStrata)
                {
                    string keyName = s1.SID.ToString() +"-"+ s2.SID.ToString();

                    if (s1.SID==s2.SID)
                        continue;
                  
                    if (s1.SPolygon.GeoPolygon.Touches(s2.SPolygon.GeoPolygon))
                    {    
                            pStrataTopoDatabase.Add(keyName, "T");
                    }
                    else
                        pStrataTopoDatabase.Add(keyName, "NT");


                }
            }
            return pStrataTopoDatabase;
        }




        /// <summary>
        /// Get all strata that contact the stratum A
        /// </summary>
        /// <param name="pStratum">stratum A</param>
        /// <param name="pTopoDatabase"></param>
        /// <returns></returns>
        public Dictionary<string, List<Stratum>> GetAllTouchStrata(Stratum pStratum, Dictionary<string, string> pTopoDatabase)
        {
            
            Dictionary<string, List<Stratum>> pTouchuStrata = new Dictionary<string, List<Stratum>>();
    
            if (pTopoDatabase == null)
                return null;

            foreach (var pS in pAllStrata)
            {
                if (pStratum.SID == pS.SID)
                    continue;
                string pStratumKey = pStratum.SID.ToString() +"-"+ pS.SID.ToString();
                
                if (pTopoDatabase[pStratumKey] == "T")
                {
                    if (pTouchuStrata.ContainsKey(pS.SName))
                    {
                        pTouchuStrata[pS.SName].Add(pS);
                    }
                    else
                    {
                        List<Stratum> pTouchStratumS = new List<Stratum>();
                        pTouchStratumS.Add(pS);
                        pTouchuStrata.Add(pS.SName, pTouchStratumS);
                    }
                }
            }

            

            return pTouchuStrata;
        }

       
    }
}
