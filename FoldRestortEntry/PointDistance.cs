using GeoAPI.Geometries;
using GeologicalEntity;
using System;
using System.Collections.Generic;
using System.Text;

namespace FoldRestortEntry
{

    /// <summary>
    /// Distance calculation between points within an R-tree
    /// </summary>
    public class PointDistance : NetTopologySuite.Index.Strtree.IItemDistance<GeoAPI.Geometries.Envelope, OccurrencePoint>
    {

        public double Distance(IBoundable<GeoAPI.Geometries.Envelope, OccurrencePoint> item1, IBoundable<GeoAPI.Geometries.Envelope, OccurrencePoint> item2)
        {
            return item1.Bounds.Distance(item2.Bounds);
        }
    }
}
