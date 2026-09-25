import os,sys
import numpy as np


"""
The following schema is for initial tests of electron vs. photon sample for
 * CNN classifier with dropout for BNN
 * Stochastic interpolant tests
 * Latent Diffusion tests
"""
SchemaDict = {}
    
def get_schema(name):
    """
    """
    from pyspark.sql.types import IntegerType, StringType, FloatType
    from petastorm.codecs import ScalarCodec, CompressedImageCodec, NdarrayCodec
    from petastorm.etl.dataset_metadata import materialize_dataset
    from petastorm.unischema import dict_to_spark_row, Unischema, UnischemaField
    from petastorm import make_reader, TransformSpec
    from petastorm.pytorch import DataLoader

    """ LOAD AND DEFINE SCHEMA HERE """
    """ We hide the pyspark and petastorm imports in case we are not using this way of storing data """
    if name not in SchemaDict:
        if name=="v0":
            SimpleDetShowerSchema_v0 = Unischema("SimpleDetShowerSchema",[
                UnischemaField('partition',np.string_, (), ScalarCodec(StringType()), False ), # label to divide the local database into folders. Allows us to split DB maker into jobs.
                UnischemaField('runid', np.int32, (), ScalarCodec(IntegerType()), False), # numerical index to divide database generation
                UnischemaField('entry', np.int64, (), ScalarCodec(IntegerType()), False), # numerical index to label entry for a given runid
                UnischemaField('pdgcode', np.int32, (), ScalarCodec(IntegerType()), False), # particle type label
                UnischemaField('depth', np.float32, (), ScalarCodec(FloatType()), False), # simulated distance from readout plane, determines transverse smearing
                UnischemaField('momentum4', np.float32, (4,), NdarrayCodec(), False), # simulated distance from readout plane, determines transverse smearing
                UnischemaField('preedeplen',np.float32, (), ScalarCodec(FloatType()), False), # distance from particle start to where first energy deposit occurs (relevant to photon)
                UnischemaField('dedx_20pix', np.float32, (20,), NdarrayCodec(), False), # dedx in first 20 pixels of start of energy deposition
                UnischemaField('edepimage', np.float32, (256,256), NdarrayCodec(), False), # energy deposition image
            ])
        else:
            raise ValueError("unrecognized schema: ",name)
    
    SchemaDict[name] = SimpleDetShowerSchema_v0
            
    return SchemaDict[name]


