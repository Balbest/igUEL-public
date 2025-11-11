from abaqus import *
from abaqusConstants import *  # for including all constants
from odbAccess import *
from odbMaterial import *
from odbSection import *
import os
#
def init_zero_tuple(input_tuple):
    # return a tuple with 0 of the same shape as input tuple
    if isinstance(input_tuple, tuple):
        return tuple(init_zero_tuple(i) for i in input_tuple)
    else:
        return 0
#
def sort_by_element(filename, output_file=None, add_zero_line=True):
    # open file
    with open(filename) as f:
        lines = f.readlines()
    elems_dic = {}
    # read line by line and find element number
    for line in lines:
        elem_number, string = line.split(">")
        elem_number = int(elem_number.replace("<", ""))
        string = string.strip()
        elems_dic[elem_number] = string

    string_final = ""
    # Sort elements by number and write to final string
    for key in sorted(elems_dic):
        string_final += " " + elems_dic[key]

    string_tuple = eval(string_final)
    string_final = str(string_tuple)[1:-1] + ","

    if add_zero_line:
        zero_tuple = init_zero_tuple(string_tuple)
        zero_tuple_str = str(zero_tuple)[1:-1] + ",\n"
        string_final = zero_tuple_str + string_final

    if output_file:
        with open(output_file, "w") as text_file:
            text_file.write(string_final)
    return string_final
#
def createODB():
#
    odb_path = "results/3D.odb"
    tmpfolder = "results/"
    #
    # ---Create an ODB (which also creates the rootAssembly)
    odb = Odb(
        name="3D_SGE_Solid_Model",
        analysisTitle="ODB created with Python ODB API",
        description="PP of IGA of 3D SGE solid",
        path=odb_path,
    )
    #
    # ---Create Material with some data, e.g., E=210000 and nu=0.3
    materialName = "Elastic Material"
    material_1 = odb.Material(name=materialName)
    material_1.Elastic(
        type=ISOTROPIC,
        temperatureDependency=OFF,
        dependencies=0,
        noCompression=OFF,
        noTension=OFF,
        moduli=LONG_TERM,
        table=((210000, 0.3),),
    )
    #
    # ---Create Section
    sectionName = "Solid Section"
    section_1 = odb.HomogeneousSolidSection(
        name=sectionName, material=materialName
    )
    #
    # ---Model data
    part1 = odb.Part(name="part-1", embeddedSpace=THREE_D, type=DEFORMABLE_BODY)
    #
    f_nodeData = open(tmpfolder + "PP_Nodes_Coords.dat", "r")
    nodeData = eval(f_nodeData.readline())
    #
    part1.addNodes(nodeData=nodeData, nodeSetName="nset-1")
    #
    f_elementData = open(tmpfolder + "PP_Elements_Nodes.dat", "r")
    elementData = eval(f_elementData.readline())
    #
    part1.addElements(elementData=elementData, type="C3D8", elementSetName="eset-1")
    #
    instance1 = odb.rootAssembly.Instance(name="part-1-1", object=part1)
    #
    # ---Create instance level sets for section assignment
    f_elLabels = open(tmpfolder + "PP_Elements.dat", "r")
    elLabels = eval(f_elLabels.readline())
    elset_1 = odb.rootAssembly.instances["part-1-1"].ElementSetFromElementLabels(
        name=materialName, elementLabels=elLabels
    )
    instance1.assignSection(region=elset_1, section=section_1)
    #
    step1 = odb.Step(
        name="step-1", description="analysis step ", domain=TIME, timePeriod=1
    )
    #
    f_Frames = open(tmpfolder + "Frames_description.dat", "r")
    f_Displacements = open(tmpfolder + "PP_U_Nodes.dat", "r")
    f_Stresses = open(tmpfolder + "PP_S_Nodes.dat", "r")
    f_Strains = open(tmpfolder + "PP_E_Nodes.dat", "r")
    f_Num_of_rows = open(tmpfolder + "PP_U_Nodes.dat", "r")
    Num = len(f_Num_of_rows.readlines())
    #
    for i, j, k, s, e in zip(
        f_Displacements, f_Frames, range(Num), f_Stresses, f_Strains
    ):
        steptime = float(j.split()[-1])
        frame = step1.Frame(incrementNumber=k, frameValue=steptime, description=j)
        #
        # -------Write nodal displacements
        uField = frame.FieldOutput(
            name="U",
            description="Displacements",
            type=VECTOR,
            validInvariants=(MAGNITUDE,)
        )
        #
        f_nodeLabelData = open(tmpfolder + "PP_Nodes.dat", "r")
        nodeLabelData = eval(f_nodeLabelData.readline())
        dispData = eval(i)
        #
        uField.addData(
            position=NODAL, instance=instance1, labels=nodeLabelData, data=dispData
        )
        #
        # -------Write nodal stresses
        sField = frame.FieldOutput(
            name="S",
            description="Stress",
            type=TENSOR_3D_FULL,
            componentLabels=("S11", "S22", "S33", "S12", "S13", "S23"),
            validInvariants=(MISES,)
        )
        #
        stressData = eval(s)
        #
        transform = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        #
        sField.addData(
            position=ELEMENT_NODAL,
            instance=instance1,
            labels=elLabels,
            data=stressData,
            localCoordSystem=transform,
        )
        #
        # -------Write nodal strains
        eField = frame.FieldOutput(
            name="E",
            description="Strain",
            type=TENSOR_3D_FULL,
            componentLabels=("E11", "E22", "E33", "E12", "E13", "E23")
        )
        #
        eData = eval(e)
        #
        transform = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        #
        eField.addData(
            position=ELEMENT_NODAL,
            instance=instance1,
            labels=elLabels,
            data=eData,
            localCoordSystem=transform,
        )
        #
        #
    if os.path.exists("results/gSx_Nodes.dat"):
       f_gx_Stresses = open(tmpfolder + "PP_gSx_Nodes.dat", "r")
       f_gy_Stresses = open(tmpfolder + "PP_gSy_Nodes.dat", "r")
       f_gz_Stresses = open(tmpfolder + "PP_gSz_Nodes.dat", "r")
       f_gx_Strains = open(tmpfolder + "PP_gEx_Nodes.dat", "r")
       f_gy_Strains = open(tmpfolder + "PP_gEy_Nodes.dat", "r")
       f_gz_Strains = open(tmpfolder + "PP_gEz_Nodes.dat", "r")
    #
       for kk, sx, sy, sz, ex, ey, ez in zip(
           range(Num), f_gx_Stresses, f_gy_Stresses, f_gz_Stresses, f_gx_Strains, f_gy_Strains, f_gz_Strains
       ):
           frame = step1.frames[kk]
           #
           # -------Write nodal double stresses, gS_xij
           gx_sField = frame.FieldOutput(
               name="gSx",
               description="Double stress (xij),",
               type=TENSOR_3D_FULL,
               componentLabels=("gSx11", "gSx22", "gSx33", "gSx12", "gSx13", "gSx23")
           )
           #
           gx_stressData = eval(sx)
           #
           transform = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
           #
           gx_sField.addData(
               position=ELEMENT_NODAL,
               instance=instance1,
               labels=elLabels,
               data=gx_stressData,
               localCoordSystem=transform,
           )
           #
           # -------Write nodal double stresses, gS_yij
           gy_sField = frame.FieldOutput(
               name="gSy",
               description="Double stress (yij),",
               type=TENSOR_3D_FULL,
               componentLabels=("gSy11", "gSy22", "gSy33", "gSy12", "gSy13", "gSy23")
           )
           #
           gy_stressData = eval(sy)
           #
           transform = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
           #
           gy_sField.addData(
               position=ELEMENT_NODAL,
               instance=instance1,
               labels=elLabels,
               data=gy_stressData,
               localCoordSystem=transform,
           )
           #
           # -------Write nodal double stresses, gS_zij
           gz_sField = frame.FieldOutput(
               name="gSz",
               description="Double stress (zij),",
               type=TENSOR_3D_FULL,
               componentLabels=("gSz11", "gSz22", "gSz33", "gSz12", "gSz13", "gSz23")
           )
           #
           gz_stressData = eval(sz)
           #
           transform = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
           #
           gz_sField.addData(
               position=ELEMENT_NODAL,
               instance=instance1,
               labels=elLabels,
               data=gz_stressData,
               localCoordSystem=transform,
           )
           #
           # -------Write nodal strain gradients, gE_xij
           gx_eField = frame.FieldOutput(
               name="gEx",
               description="Strain gradient (xij),",
               type=TENSOR_3D_FULL,
               componentLabels=("gEx11", "gEx22", "gEx33", "gEx12", "gEx13", "gEx23")
           )
           #
           gx_strainData = eval(ex)
           #
           transform = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
           #
           gx_eField.addData(
               position=ELEMENT_NODAL,
               instance=instance1,
               labels=elLabels,
               data=gx_strainData,
               localCoordSystem=transform,
           )
           #
           # -------Write nodal strain gradients, gE_yij
           gy_eField = frame.FieldOutput(
               name="gEy",
               description="Strain gradient (yij),",
               type=TENSOR_3D_FULL,
               componentLabels=("gEy11", "gEy22", "gEy33", "gEy12", "gEy13", "gEy23")
           )
           #
           gy_strainData = eval(ey)
           #
           transform = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
           #
           gy_eField.addData(
               position=ELEMENT_NODAL,
               instance=instance1,
               labels=elLabels,
               data=gy_strainData,
               localCoordSystem=transform,
           )
           #
           # -------Write nodal strain gradients, gE_zij
           gz_eField = frame.FieldOutput(
               name="gEz",
               description="Strain gradient (zij),",
               type=TENSOR_3D_FULL,
               componentLabels=("gEz11", "gEz22", "gEz33", "gEz12", "gEz13", "gEz23")
           )
           #
           gz_strainData = eval(ez)
           #
           transform = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
           #
           gz_eField.addData(
               position=ELEMENT_NODAL,
               instance=instance1,
               labels=elLabels,
               data=gz_strainData,
               localCoordSystem=transform,
           )
           #
# -----------------------------------------------------------------------------------------
#
    odb.save()
    odb.close()
#
if __name__ == "__main__":
    #
    # ###################################################################
    # Assembling Element, Node, Displacement, Stress and Strain data
    # ###################################################################
    #
    files = [
        "Nodes.dat",
        "Nodes_Coords.dat",
        "Elements_Nodes.dat",
        "Elements.dat",
    ]
#
    for file in files:
        sort_by_element("results/" + file, "results/PP_" + file, False)
#
    files = [
        "E_Nodes.dat",
        "S_Nodes.dat",
        "U_Nodes.dat",
    ]
#
    for file in files:
        sort_by_element("results/" + file, "results/PP_" + file, True)
#
    if os.path.exists("results/gSx_Nodes.dat"):
       files = [
           "gSx_Nodes.dat",
           "gSy_Nodes.dat",
           "gSz_Nodes.dat",
           "gEx_Nodes.dat",
           "gEy_Nodes.dat",
           "gEz_Nodes.dat",
       ]
#
       for file in files:
           sort_by_element("results/" + file, "results/PP_" + file, True)
#
    createODB()