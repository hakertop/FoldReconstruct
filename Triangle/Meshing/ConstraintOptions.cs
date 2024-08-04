
namespace TriangleNet.Meshing
{
    /// <summary>
    /// Mesh constraint options for polygon triangulation.
    /// </summary>
    public class ConstraintOptions
    {
        // TODO: remove ConstraintOptions.UseRegions

        /// <summary>
        /// Gets or sets a value indicating whether to use regions.
        /// 获取或设置一个值，该值指示是否使用区域。
        /// </summary>
        [System.Obsolete("Not used anywhere, will be removed in beta 4.")]
        public bool UseRegions { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether to create a Conforming
        /// Delaunay triangulation.
        /// 获取或设置一个值，该值指示是否创建符合Delaunay三角剖分
        /// </summary>
        public bool ConformingDelaunay { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether to enclose the convex
        /// hull with segments.
        /// 获取或设置一个值，该值指示是否用段包围凸包。
        /// </summary>
        public bool Convex { get; set; }

        /// <summary>
        /// Gets or sets a flag indicating whether to suppress boundary
        /// segment splitting.
        /// 获取或设置一个标志，该标志指示是否抑制边界段分割。
        /// </summary>
        /// <remarks>
        /// 0 = split segments (default)
        /// 0 = 分离段
        /// 1 = no new vertices on the boundary
        /// 1 = 边界上没有新的顶点
        /// 2 = prevent all segment splitting, including internal boundaries
        /// 2 = 防止所有的段分裂，包括内部边界
        /// </remarks>
        public int SegmentSplitting { get; set; }
    }
}
