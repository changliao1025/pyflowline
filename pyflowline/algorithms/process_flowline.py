import os, sys

# will use logger for better logging in the future
from pyflowline.formats.export_flowline import export_flowline_to_geojson
from pyflowline.formats.export_vertex import export_vertex_to_geojson
from pyflowline.algorithms.split.split_by_length import split_flowline_by_length
from pyflowline.classes.timer import pytimer
from pyearthriver.classes.rivergraph import pyrivergraph


def process_flowline(
    aFlowline_in,
    pVertex_outlet_initial,
    sWorkspace_output_basin,
    iFlag_remove_duplicate_river=0,
    iFlag_remove_disconnected_river=0,
    iFlag_merge_flowline=1,
    iFlag_remove_braided_river=0,
    iFlag_remove_parallel_river=0,
    iFlag_remove_cycle=0,
    iFlag_remove_small_river=0,
    iFlag_break_by_distance=0,
    iFlag_direction_insensitive=0,
    iFlag_debug=0,
    dThreshold_small_river=10000,
    dThreshold_break_by_distance=5000.0,
):

    ptimer = pytimer()
    try:
        print("Find flowline vertex")
        ptimer.start()
        pRivergraph = pyrivergraph(aFlowline_in, pVertex_outlet_initial)
        aVertex = pRivergraph.aVertex
        ptimer.stop()
        if iFlag_debug == 1:
            sFilename_out = (
                "flowline_vertex_without_confluence_before_intersect.geojson"
            )
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_vertex_to_geojson(aVertex, sFilename_out)
    except:
        print("Error in find flowline vertex")

    sys.stdout.flush()

    try:
        print("Start splitting flowline")
        ptimer.start()
        aFlowline_basin_simplified = pRivergraph.split_flowline()
        ptimer.stop()
        if iFlag_debug == 1:
            sTime = ptimer.now()
            sFilename_out = "flowline_split_by_point_" + sTime + ".geojson"
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_geojson(aFlowline_basin_simplified, sFilename_out)
    except:
        print("Error in split flowline")

    sys.stdout.flush()

    if iFlag_remove_duplicate_river == 1:
        try:
            print("Start removing duplicate flowline")
            ptimer.start()
            aFlowline_basin_conceptual = pRivergraph.remove_duplicate_flowlines(
                iFlag_direction_insensitive=iFlag_direction_insensitive
            )
            ptimer.stop()
            if iFlag_debug == 1:
                sTime = ptimer.now()
                sFilename_out = "flowline_remove_duplicate_" + sTime + ".geojson"
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_basin_conceptual, sFilename_out)
        except:
            print("Error in remove_duplicate_flowline.")

    sys.stdout.flush()

    if iFlag_remove_disconnected_river == 1:
        try:
            print("Start removing disconnected flowline")
            ptimer.start()
            aFlowline_basin_conceptual = pRivergraph.remove_disconnected_flowlines()
            ptimer.stop()
            if iFlag_debug == 1:
                sTime = ptimer.now()
                sFilename_out = "flowline_remove_disconnected_" + sTime + ".geojson"
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_basin_conceptual, sFilename_out)
        except:
            print("Error in remove disconnected flowlines.")

    sys.stdout.flush()

    if iFlag_merge_flowline == 1:
        try:
            print("Start merging flowline")
            ptimer.start()
            aFlowline_basin_conceptual = pRivergraph.merge_flowline()
            ptimer.stop()
            if iFlag_debug == 1:
                sTime = ptimer.now()
                sFilename_out = "flowline_merge_flowline_" + sTime + ".geojson"
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_basin_conceptual, sFilename_out)
        except:
            print("Error in merge_flowline.")

    if iFlag_remove_braided_river == 1:
        try:
            print("Start removing braided river")
            ptimer.start()
            aFlowline_basin_simplified = pRivergraph.remove_braided_river()
            ptimer.stop()
            if iFlag_debug == 1:
                sTime = ptimer.now()
                sFilename_out = "flowline_remove_braided_river_" + sTime + ".geojson"
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_basin_simplified, sFilename_out)
        except:
            print("Error in brainded river removal")

    sys.stdout.flush()

    if iFlag_remove_parallel_river == 1:
        try:
            print("Start removing parallel river")
            ptimer.start()
            aFlowline_basin_simplified = pRivergraph.remove_parallel_river()
            ptimer.stop()
            if iFlag_debug == 1:
                sTime = ptimer.now()
                sFilename_out = "flowline_remove_parallel_" + sTime + ".geojson"
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_basin_simplified, sFilename_out)
        except:
            print("Error in parallel river removal")
        sys.stdout.flush()

    if iFlag_remove_cycle == 1:
        try:
            print("Start cycle removal")
            ptimer.start()
            aFlowline_basin_conceptual = pRivergraph.remove_cycle()
            ptimer.stop()
            if iFlag_debug == 1:
                sTime = ptimer.now()
                sFilename_out = "flowline_remove_loop_" + sTime + ".geojson"
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_basin_conceptual, sFilename_out)
        except:
            print("Error in cycle removal.")

    try:
        print("Started update stream order initial")
        ptimer.start()
        aFlowline_basin_simplified = pRivergraph.update_headwater_stream_order()
        pRivergraph.define_river_confluence()
        pRivergraph.define_stream_order(iFlag_so_method_in=1)
        ptimer.stop()
    except:
        print("Error in update stream order initial")

    sys.stdout.flush()

    if iFlag_remove_small_river == 1:
        try:
            print("Started small river removal")
            ptimer.start()
            aFlowline_basin_simplified = pRivergraph.remove_small_river(
                dThreshold_small_river=dThreshold_small_river,
                nIterations=3,
                iFlag_debug=iFlag_debug,
                sWorkspace_output_basin=sWorkspace_output_basin,
            )
            ptimer.stop()
        except:
            print("Error in small river removal")
    else:
        try:
            print("Basin ", "started flowline merging (no small river removal)")
            ptimer.start()
            aFlowline_basin_simplified = pRivergraph.merge_flowline()
            aFlowline_basin_simplified = pRivergraph.update_headwater_stream_order()
            ptimer.stop()
            if iFlag_debug == 1:
                sFilename_out = "flowline_merge_before_intersect.geojson"
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(aFlowline_basin_simplified, sFilename_out)
        except:
            print("Error in flowline merging (no small river removal)")

    sys.stdout.flush()

    try:
        print("Started stream segment definition")
        ptimer.start()
        aFlowline_basin_simplified, aStream_segment = (
            pRivergraph.define_stream_segment()
        )
        ptimer.stop()
        if iFlag_debug == 1:
            sTime = ptimer.now()
            sFilename_out = os.path.join(
                sWorkspace_output_basin, "flowline_segment_index_" + sTime + ".geojson"
            )
            export_flowline_to_geojson(
                aFlowline_basin_simplified,
                sFilename_out,
                aAttribute_data=[aStream_segment],
                aAttribute_field=["stream_segment"],
                aAttribute_dtype=["int"],
            )
    except:
        print("Error in stream segment definition")
    sys.stdout.flush()

    try:
        print("Started confluence definition")
        ptimer.start()
        aConfluence_basin_simplified = pRivergraph.define_river_confluence()
        ptimer.stop()
    except:
        print("Error in confluence definition")
    sys.stdout.flush()

    try:
        print("Started stream topology definition")
        ptimer.start()
        aFlowline_basin_simplified = pRivergraph.define_stream_topology()
        ptimer.stop()
    except:
        print("Error in stream topology definition")
    sys.stdout.flush()

    try:
        print("Started stream order definition")
        ptimer.start()
        aFlowline_basin_simplified, aStream_order = pRivergraph.define_stream_order()
        ptimer.stop()
    except:
        print("Error in stream order definition")

    sys.stdout.flush()

    if iFlag_break_by_distance == 1:
        try:
            print("Started flowline split by length")
            ptimer.start()
            aFlowline_basin_simplified_split = split_flowline_by_length(
                aFlowline_basin_simplified, dThreshold_break_by_distance
            )
            ptimer.stop()
            if iFlag_debug == 1:
                sTime = ptimer.now()
                sFilename_out = "flowline_split_by_length_" + sTime + ".geojson"
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_geojson(
                    aFlowline_basin_simplified_split, sFilename_out
                )
        except:
            print("Error in flowline split by length")

    return (
        aFlowline_basin_simplified,
        aConfluence_basin_simplified,
        aStream_segment,
        aStream_order,
    )
