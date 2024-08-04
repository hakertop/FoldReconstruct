using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GeologicalEntity
{
    //define the fold struct
    public class Fold
    {
        //field

        public List<int> foldpolindex = new List<int>();//All stratigraphic sequences traversed by the fold, with numbers corresponding to the geological units' sequence in the geological map
        public List<Stratum> foldpolygon = new List<Stratum>();//All stratigraphic layers encompassed by the fold
        public int coreid;//Sequence number of the core stratigraphic layer
        public FoldType foldtype;//fold type
        public SumFoldType sumfoldtype;//syncline or anticline
        public double length;//the long axis of the fold
        public double width; //the short axis of the fold
        public Hinge hingle;//the attitude of the hingle
        public AxialPlane axialplane;//the attitude of the axial plane
        public LeftLimb leftloccurrence;//dominant attitude on the left side
        public RightLimb rightloccurrence;//dominant attitude on the right side
        public int limbangle;//interlimb
        public TypeOfLimbAngle typeoflimbangle;//interlimb type
        public int sectionlinenum;//index of the cross-section line
    }

    /// <summary>
    /// fold types，can be classfied into seven types
    /// </summary>
    public enum FoldType
    {
        VerticalHorizontalFold,
        VerticalPlungingFold,
        PourverticaloldF,
        ReclinedFold,
        InclinedHorizontalFold,
        RecumbentFold,
        InclinedPlungingFold
    }

    /// <summary>
    /// limbangle type
    /// </summary>
    public enum TypeOfLimbAngle
    {
        Gently,
        Open,
        Close,
        Tight,
        Inclined
    }

    public enum SumFoldType
    {
        Syncline,
        Anticline,
        NotFold
    }


    /// <summary>
    /// the attitude of the hinge
    /// </summary>
    public class Hinge
    {
        public int HdipAngle;//dip angle
        public int Htendency;//dip
    }

    /// <summary>
    /// the attitude of the axial plane
    /// </summary>
    public class AxialPlane
    {
        public int AdipAngle;//dip angle
        public int Atendency;//dip
    }

    /// <summary>
    /// the dominant attitude of the left limb
    /// </summary>
    public class LeftLimb
    {
        public int LdipAngle;//dip angle
        public int Ltendency;//dip
    }

    /// <summary>
    /// the dominant attitude of the right limb
    /// </summary>
    public class RightLimb
    {
        public int RdipAngle;//dip angle
        public int Rtendency;//dip
    }

}
