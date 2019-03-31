from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (QgsField, QgsFields, QgsFeature, QgsGeometry, QgsFeatureSink, QgsFeatureRequest, QgsProcessing, QgsProcessingAlgorithm, QgsProcessingParameterFeatureSource, QgsProcessingOutputNumber, QgsProcessingParameterField, QgsProcessingParameterBoolean, QgsProcessingException)
                       
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


    class SequenceGenerator():
        def __init__(self,centroid_layer,trajectory_layer,weight_field=None):
            centroids = [f for f in centroid_layer.getFeatures()]
            self.cell_index = QgsSpatialIndex()
            for f in centroids:
                self.cell_index.insertFeature(f)
            self.id_to_centroid = {f.id(): [f,[0,0,0,0,0]] for (f) in centroids}
            self.weight_field = weight_field
            if weight_field is not None:
                self.weightIdx = trajectory_layer.fieldNameIndex(weight_field)
            else:
                self.weightIdx = None
            self.sequences = {}
            
            nTraj = float(trajectory_layer.featureCount())
            for i,traj in enumerate(trajectory_layer.getFeatures()):
                self.evaluate_trajectory(traj)
                progress.setPercentage(i/nTraj*100)
                
        def evaluate_trajectory(self,trajectory):
            points = trajectory.geometry().asPolyline()
            this_sequence = []
            for i, pt in enumerate(points):
                id = self.cell_index.nearestNeighbor(pt,1)[0]
                nearest_cell = self.id_to_centroid[id][0]
                nearest_cell_id = nearest_cell.id()
                prev_cell_id = None
                if len(this_sequence) > 1:
                    prev_cell_id = this_sequence[-1]
                    if self.weight_field is not None:
                        weight = trajectory.attributes()[self.weightIdx]
                    else:
                        weight = 1
                    if self.sequences.has_key((prev_cell_id,nearest_cell_id)):
                        self.sequences[(prev_cell_id,nearest_cell_id)] += weight
                    else:
                        self.sequences[(prev_cell_id,nearest_cell_id)] = weight
                if nearest_cell_id != prev_cell_id: 
                    # we have changed to a new cell --> up the counter 
                    m = trajectory.geometry().geometry().pointN(i).m()
                    t = datetime(1970,1,1) + timedelta(seconds=m) + timedelta(hours=8) # Beijing GMT+8
                    h = t.hour 
                    self.id_to_centroid[id][1][0] = self.id_to_centroid[id][1][0] + 1
                    self.id_to_centroid[id][1][h/6+1] = self.id_to_centroid[id][1][h/6+1] + 1
                    this_sequence.append(nearest_cell_id)
        
        def create_flow_lines(self):
            lines = []
            for key,value in self.sequences.iteritems(): 
                p1 = self.id_to_centroid[key[0]][0].geometry().asPoint()
                p2 = self.id_to_centroid[key[1]][0].geometry().asPoint()
                feat = QgsFeature()
                feat.setGeometry(QgsGeometry.fromPolyline([p1,p2]))
                feat.setAttributes([key[0],key[1],value])
                lines.append(feat)
            return lines


    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_TRAJECTORIES,
            self.tr(self.INPUT_TRAJECTORIES),
            [QgsProcessing.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterField(self.WEIGHT_FIELD,
            self.tr(self.WEIGHT_FIELD),
            'Weight',
            self.INPUT_TRAJECTORIES,
            QgsProcessingParameterField.Numeric))
        self.addParameter(QgsProcessingParameterBoolean(self.USE_WEIGHT_FIELD,
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
        flowIdx = trajectory_layer.fields().indexFromName(weight_field)
        use_weight_field = self.parameterAsBool(parameters, self.USE_WEIGHT_FIELD, context)
        lineFields = QgsFields()
        lineFields.append(QgsField('FROM', QVariant.Int))
        lineFields.append(QgsField('TO', QVariant.Int))
        lineFields.append(QgsField('COUNT', QVariant.Int))
        (lineSink, line_dest_id) = self.parameterAsSink(parameters, self.OUTPUT_FLOWLINES, context,
                                        lineFields, trajectory_layer.wkbType(), trajectory_layer.sourceCrs())        
        pointFields = centroid_layer.fields()
        pointFields.append(QgsField('COUNT',QVariant.Int))
        pointFields.append(QgsField('COUNT_Q1',QVariant.Int))
        pointFields.append(QgsField('COUNT_Q2',QVariant.Int))
        pointFields.append(QgsField('COUNT_Q3',QVariant.Int))
        pointFields.append(QgsField('COUNT_Q4',QVariant.Int))
        (pointSink, point_dest_id) = self.parameterAsSink(parameters, self.OUTPUT_CELL_COUNTS, context,
                                        pointFields, centroid_layer.wkbType(), centroid_layer.sourceCrs()) 
        
        sg = SequenceGenerator(centroid_layer,trajectory_layer, weight_field if use_weight_field else None)
        
        geom_type = 2
        for f in sg.create_flow_lines():
            lineSink.addFeature(f, QgsFeatureSink.FastInsert)
        
        for key, value in sg.id_to_centroid.iteritems():
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

