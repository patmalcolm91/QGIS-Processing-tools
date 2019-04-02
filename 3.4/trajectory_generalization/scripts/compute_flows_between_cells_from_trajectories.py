from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (QgsSpatialIndex, QgsWkbTypes, QgsField, QgsFields, QgsFeature, QgsGeometry, QgsPoint,
                       QgsFeatureSink, QgsProcessingParameterFeatureSink, QgsProcessing, QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource, QgsProcessingParameterField, QgsProcessingParameterBoolean,
                       QgsProcessingParameterNumber)
from datetime import datetime, timedelta
import math
import time

               
class ComputeFlowsBetweenCellsFromTrajectories(QgsProcessingAlgorithm):
    INPUT_TRAJECTORIES = 'Input Trajectories'
    WEIGHT_FIELD = 'Weight Field'
    USE_WEIGHT_FIELD = 'Use weight field'
    INPUT_CELL_CENTERS = 'Input Cell Centers'
    OUTPUT_FLOWLINES = 'Output Flowlines'
    OUTPUT_CELL_COUNTS = 'Output Cell Counts'
    TIMEZONE = "Timezone (UTC+...)"
 
    def __init__(self):
        super().__init__()
 
    def name(self):
        return "Compute Flows Between Cells fom Trajectories"
     
    def tr(self, text):
        return QCoreApplication.translate("ComputeFlowsBetweenCellsFromTrajectories", text)
         
    def displayName(self):
        return self.tr("Compute Flows Between Cells fom Trajectories")
 
    def group(self):
        return self.tr("Trajectory Generalization")
 
    def groupId(self):
        return "trajectoryGeneralization"
 
    def shortHelpString(self):
        return self.tr("""
        Computes flows between cells from trajectories.
        Parameters:
        &bull; Input Trajectories: the LineM or Line layer containing the trajectories to be aggregated.
        Note: if using OD Lines as input, be sure to first densify the lines at an appropriate interval.
        &bull; Weight Field: Field in the trajectory layer which contains the weight of each trajectory. Default False.
        &bull; Use Weight Field: whether or not to use the specified weight field. If False, a weight of 1 is used.
        &bull; Input Cell Centers: the points to which the trajectories should be snapped before calculating.
        &bull; Timezone: offset of the data from UTC (for splitting counts in to 6 hour windows). Default local offset.
        &bull; Output Flowlines: line layer containing the total flows between each connected node.
        &bull; Output Cell Counts: point layer containing the total flow through each point, also divided in 6h windows.
        """)
 
    def helpUrl(self):
        return "https://qgis.org"
         
    def createInstance(self):
        return type(self)()

    class SequenceGenerator:
        def __init__(self, centroid_layer, trajectory_layer, feedback, timezone, weight_field=None):
            centroids = [f for f in centroid_layer.getFeatures()]
            self.cell_index = QgsSpatialIndex()
            for f in centroids:
                self.cell_index.insertFeature(f)
            self.id_to_centroid = {f.id(): [f, [0, 0, 0, 0, 0]] for (f) in centroids}
            self.timezone = timezone
            self.weight_field = weight_field
            if weight_field is not None:
                self.weightIdx = trajectory_layer.fields().indexFromName(weight_field)
            else:
                self.weightIdx = None
            self.sequences = {}
            
            n_traj = float(trajectory_layer.featureCount())
            for i, traj in enumerate(trajectory_layer.getFeatures()):
                self.evaluate_trajectory(traj)
                feedback.setProgress(i/n_traj*100)
                
        def evaluate_trajectory(self, trajectory):
            points = trajectory.geometry().asPolyline()
            this_sequence = []
            weight = 1 if self.weight_field is None else trajectory.attributes()[self.weightIdx]
            prev_cell_id = None
            for i, pt in enumerate(points):
                nn_id = self.cell_index.nearestNeighbor(pt, 1)[0]
                nearest_cell = self.id_to_centroid[nn_id][0]
                nearest_cell_id = nearest_cell.id()
                if len(this_sequence) >= 1:
                    prev_cell_id = this_sequence[-1]
                    if nearest_cell_id != prev_cell_id:
                        if (prev_cell_id, nearest_cell_id) in self.sequences:
                            self.sequences[(prev_cell_id, nearest_cell_id)] += weight
                        else:
                            self.sequences[(prev_cell_id, nearest_cell_id)] = weight
                if nearest_cell_id != prev_cell_id: 
                    # we have changed to a new cell --> up the counter 
                    m = trajectory.geometry().vertexAt(i).m()
                    if math.isnan(m):
                        m = 0
                    t = datetime(1970, 1, 1) + timedelta(seconds=m) + timedelta(hours=self.timezone)
                    h = t.hour
                    self.id_to_centroid[nn_id][1][0] += weight
                    self.id_to_centroid[nn_id][1][int(h/6)+1] += weight
                    this_sequence.append(nearest_cell_id)
        
        def create_flow_lines(self):
            lines = []
            for key, value in self.sequences.items():
                p1 = self.id_to_centroid[key[0]][0].geometry().asPoint()
                p2 = self.id_to_centroid[key[1]][0].geometry().asPoint()
                p1 = QgsPoint(p1.x(), p1.y())
                p2 = QgsPoint(p2.x(), p2.y())
                feat = QgsFeature()
                feat.setGeometry(QgsGeometry.fromPolyline([p1, p2]))
                feat.setAttributes([key[0], key[1], value])
                lines.append(feat)
            return lines

    def initAlgorithm(self, config=None):
        local_timezone = time.timezone if (time.localtime().tm_isdst == 0) else time.altzone
        local_timezone /= -3600  # gets the time zone offset of the local machine to use a default

        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_TRAJECTORIES,
            self.tr(self.INPUT_TRAJECTORIES),
            [QgsProcessing.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterField(
            self.WEIGHT_FIELD,
            self.tr(self.WEIGHT_FIELD),
            'Weight',
            self.INPUT_TRAJECTORIES,
            QgsProcessingParameterField.Numeric))
        self.addParameter(QgsProcessingParameterBoolean(
            self.USE_WEIGHT_FIELD,
            self.tr(self.USE_WEIGHT_FIELD),
            QVariant(False)))
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_CELL_CENTERS,
            self.tr(self.INPUT_CELL_CENTERS),
            [QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_FLOWLINES,
            self.tr(self.OUTPUT_FLOWLINES),
            QgsProcessing.TypeVectorLine))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_CELL_COUNTS,
            self.tr(self.OUTPUT_CELL_COUNTS),
            QgsProcessing.TypeVectorPoint))
        self.addParameter(QgsProcessingParameterNumber(
            self.TIMEZONE,
            self.tr(self.TIMEZONE),
            QgsProcessingParameterNumber.Integer,
            QVariant(local_timezone),
            False,
            -12,
            12
        ))

    def processAlgorithm(self, parameters, context, feedback):
        timezone = self.parameterAsInt(parameters, self.TIMEZONE, context)
        centroid_layer = self.parameterAsSource(parameters, self.INPUT_CELL_CENTERS, context)
        trajectory_layer = self.parameterAsSource(parameters, self.INPUT_TRAJECTORIES, context)
        weight_field = self.parameterAsString(parameters, self.WEIGHT_FIELD, context)
        use_weight_field = self.parameterAsBool(parameters, self.USE_WEIGHT_FIELD, context)
        line_fields = QgsFields()
        line_fields.append(QgsField('FROM', QVariant.Int))
        line_fields.append(QgsField('TO', QVariant.Int))
        line_fields.append(QgsField('COUNT', QVariant.Int))
        (lineSink, line_dest_id) = self.parameterAsSink(parameters, self.OUTPUT_FLOWLINES, context, line_fields,
                                                        QgsWkbTypes.LineString, trajectory_layer.sourceCrs())
        point_fields = centroid_layer.fields()
        point_fields.append(QgsField('COUNT', QVariant.Int))
        point_fields.append(QgsField('COUNT_Q1', QVariant.Int))
        point_fields.append(QgsField('COUNT_Q2', QVariant.Int))
        point_fields.append(QgsField('COUNT_Q3', QVariant.Int))
        point_fields.append(QgsField('COUNT_Q4', QVariant.Int))
        (pointSink, point_dest_id) = self.parameterAsSink(parameters, self.OUTPUT_CELL_COUNTS, context, point_fields,
                                                          centroid_layer.wkbType(), centroid_layer.sourceCrs())
        
        sg = self.SequenceGenerator(centroid_layer, trajectory_layer, feedback, timezone,
                                    weight_field if use_weight_field else None)

        for f in sg.create_flow_lines():
            lineSink.addFeature(f, QgsFeatureSink.FastInsert)
        
        for key, value in sg.id_to_centroid.items():
            (in_feature, n) = value
            out_feature = QgsFeature()
            out_feature.setGeometry(in_feature.geometry())
            attributes = in_feature.attributes()
            attributes.append(n[0])
            attributes.append(n[1])
            attributes.append(n[2])
            attributes.append(n[3])
            attributes.append(n[4])
            out_feature.setAttributes(attributes)
            pointSink.addFeature(out_feature, QgsFeatureSink.FastInsert)
 
        return {self.OUTPUT_FLOWLINES: line_dest_id,
                self.OUTPUT_CELL_COUNTS: point_dest_id}
