from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (QgsSpatialIndex, QgsWkbTypes, QgsField, QgsFields, QgsFeature, QgsGeometry, QgsPoint,
                       QgsFeatureSink, QgsProcessingParameterFeatureSink, QgsProcessing, QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource, QgsProcessingParameterField, QgsProcessingParameterBoolean)
from datetime import datetime, timedelta
import math

               
class ComputeFlowsBetweenCellsFromTrajectories(QgsProcessingAlgorithm):
    INPUT_TRAJECTORIES = 'Input Trajectories'
    WEIGHT_FIELD = 'Weight Field'
    USE_WEIGHT_FIELD = 'Use weight field'
    INPUT_CELL_CENTERS = 'Input Cell Centers'
    OUTPUT_FLOWLINES = 'Output Flowlines'
    OUTPUT_CELL_COUNTS = 'Output Cell Counts'
 
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
        """)
 
    def helpUrl(self):
        return "https://qgis.org"
         
    def createInstance(self):
        return type(self)()

    class SequenceGenerator:
        def __init__(self, centroid_layer, trajectory_layer, feedback, weight_field=None):
            centroids = [f for f in centroid_layer.getFeatures()]
            self.cell_index = QgsSpatialIndex()
            for f in centroids:
                self.cell_index.insertFeature(f)
            self.id_to_centroid = {f.id(): [f, [0, 0, 0, 0, 0]] for (f) in centroids}
            self.weight_field = weight_field
            if weight_field is not None:
                self.weightIdx = trajectory_layer.fields().indexFromName(weight_field)
            else:
                self.weightIdx = None
            self.sequences = {}
            
            nTraj = float(trajectory_layer.featureCount())
            for i, traj in enumerate(trajectory_layer.getFeatures()):
                self.evaluate_trajectory(traj)
                feedback.setProgress(i/nTraj*100)
                
        def evaluate_trajectory(self,trajectory):
            points = trajectory.geometry().asPolyline()
            this_sequence = []
            weight = 1 if self.weight_field is None else trajectory.attributes()[self.weightIdx]
            prev_cell_id = None
            for i, pt in enumerate(points):
                id = self.cell_index.nearestNeighbor(pt,1)[0]
                nearest_cell = self.id_to_centroid[id][0]
                nearest_cell_id = nearest_cell.id()
                if len(this_sequence) >= 1:
                    prev_cell_id = this_sequence[-1]
                    if nearest_cell_id != prev_cell_id:
                        if (prev_cell_id,nearest_cell_id) in self.sequences:
                            self.sequences[(prev_cell_id,nearest_cell_id)] += weight
                        else:
                            self.sequences[(prev_cell_id,nearest_cell_id)] = weight
                if nearest_cell_id != prev_cell_id: 
                    # we have changed to a new cell --> up the counter 
                    m = trajectory.geometry().vertexAt(i).m()
                    if math.isnan(m):
                        m = 0
                    t = datetime(1970, 1, 1) + timedelta(seconds=m) + timedelta(hours=8)  # Beijing GMT+8
                    h = t.hour
                    self.id_to_centroid[id][1][0] += weight
                    self.id_to_centroid[id][1][int(h/6)+1] += weight
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

    def processAlgorithm(self, parameters, context, feedback):
        centroid_layer = self.parameterAsSource(parameters, self.INPUT_CELL_CENTERS, context)
        trajectory_layer = self.parameterAsSource(parameters, self.INPUT_TRAJECTORIES, context)
        weight_field = self.parameterAsString(parameters, self.WEIGHT_FIELD, context)
        use_weight_field = self.parameterAsBool(parameters, self.USE_WEIGHT_FIELD, context)
        lineFields = QgsFields()
        lineFields.append(QgsField('FROM', QVariant.Int))
        lineFields.append(QgsField('TO', QVariant.Int))
        lineFields.append(QgsField('COUNT', QVariant.Int))
        (lineSink, line_dest_id) = self.parameterAsSink(parameters, self.OUTPUT_FLOWLINES, context,
                                        lineFields, QgsWkbTypes.LineString, trajectory_layer.sourceCrs())
        pointFields = centroid_layer.fields()
        pointFields.append(QgsField('COUNT',QVariant.Int))
        pointFields.append(QgsField('COUNT_Q1',QVariant.Int))
        pointFields.append(QgsField('COUNT_Q2',QVariant.Int))
        pointFields.append(QgsField('COUNT_Q3',QVariant.Int))
        pointFields.append(QgsField('COUNT_Q4',QVariant.Int))
        (pointSink, point_dest_id) = self.parameterAsSink(parameters, self.OUTPUT_CELL_COUNTS, context,
                                        pointFields, centroid_layer.wkbType(), centroid_layer.sourceCrs()) 
        
        sg = self.SequenceGenerator(centroid_layer, trajectory_layer, feedback,
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

