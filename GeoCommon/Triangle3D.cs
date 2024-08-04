using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GeoCommon
{
    public class Triangle3D
    {

        /// <summary>
        /// trangle ID
        /// </summary>
        public int id;

        /// <summary>
        /// trangle 0 ID
        /// </summary>
        public int v0;

        /// <summary>
        /// trangle 1 ID
        /// </summary>
        public int v1;

        /// <summary>
        /// trangle 2 ID
        /// </summary>
        public int v2;

        /// <summary>
        /// 
        /// </summary>
        public int triMeshId;

        /// <summary>
        /// 
        /// </summary>
        public TriangleType type = TriangleType.Unknown;

        /// <summary>
        /// 
        /// </summary>
        public Triangle3D()
        {
        }
    }

    /// <summary>
    /// 
    /// </summary>
    public enum TriangleType
    {
        Unknown,
        Normal,
        Terrian,
        Bundary
    }
}
