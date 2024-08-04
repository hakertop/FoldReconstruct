
namespace TriangleNet.Meshing
{
    using System;
    using TriangleNet.Geometry;

    /// <summary>
    /// Mesh constraint options for quality triangulation.
    /// </summary>
    public class QualityOptions
    {
        /// <summary>
        /// Gets or sets a maximum angle constraint.
        /// 获取或设置最大角度约束。
        /// </summary>
        public double MaximumAngle { get; set; }

        /// <summary>
        /// Gets or sets a minimum angle constraint.
        /// 获取或设置最小角度约束。
        /// </summary>
        public double MinimumAngle { get; set; }

        /// <summary>
        /// Gets or sets a maximum triangle area constraint.
        /// 获取或设置最大三角形面积约束。
        /// </summary>
        public double MaximumArea { get; set; }

        /// <summary>
        /// 获取或设置用户定义的三角形约束。
        /// </summary>
        /// <remarks>
        /// The test function will be called for each triangle in the mesh. The
        /// second argument is the area of the triangle tested. If the function
        /// returns true, the triangle is considered bad and will be refined.
        /// 测试函数将为网格中的每个三角形调用。第二个参数是测试的三角形的面积。
        /// 如果函数返回true，则认为三角形是坏的，将被细化。
        /// </remarks>
        public Func<ITriangle, double, bool> UserTest { get; set; }

        /// <summary>
        /// Gets or sets an area constraint per triangle.
        /// 获取或设置每个三角形的面积约束。
        /// </summary>
        /// <remarks>
        /// If this flag is set to true, the <see cref="ITriangle.Area"/> value will
        /// be used to check if a triangle needs refinement.
        /// 如果此标志设置为true，则<see cref="三角。Area"/>值将用于检查一个三角形是否需要细化。
        /// </remarks>
        public bool VariableArea { get; set; }

        /// <summary>
        /// Gets or sets the maximum number of Steiner points to be inserted into the mesh.
        /// 获取或设置要插入网格中的斯坦纳点的最大数目。
        /// </summary>
        /// <remarks>
        /// If the value is 0 (default), an unknown number of Steiner points may be inserted
        /// to meet the other quality constraints.
        /// 如果该值为0(默认值)，则可以插入未知数量的斯坦纳点，以满足其他质量约束。
        /// </remarks>
        public int SteinerPoints { get; set; }
    }
}
