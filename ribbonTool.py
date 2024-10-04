import maya.cmds as cmds
from builtins import zip

model_tool = None
first_tool = None
addFk = None
addDeformation = None

def select_model(*args):
    select = cmds.ls(selection=True, flatten=True)
    model_str = ",".join(select)
    cmds.textField(model_tool, edit = True, text = model_str)

def select_edges(*args):
    selected_edges = cmds.ls(selection = True, flatten = True)
    edges_str = ",".join(selected_edges)
    cmds.textField(first_tool, edit = True, text = edges_str)

def createFol(num_edges, loft, select_geo):
    fol_grp = cmds.group(em = True, name = "follicle_grp")
    pV_num = 0.0

    ribbon_transform_grp = cmds.group(em = True, name = "ribbon_transform_grp")
    ribbon_noTransform_grp = cmds.group(em = True, n = "ribbon_noTransform_grp")
    cmds.hide(loft)
    cmds.parent(loft, ribbon_transform_grp)

    joints = []
    ctrl_grp = []

    created_circle_root = {}
    created_joint_root = {}

    ctrl_dict = {}
    ctrl_constraint_dict = {}

    for num in range(num_edges + 1):
        follicle_name = "follicle_{0}".format(num)
        fol_ctrl_grp_name = "fol_ctrl_grp_{0}".format(num)
        fol_ctrl_name = "fol_ctrl_{0}".format(num)
        joint_name = "joint_{0}".format(num)
        
        follicleShp = cmds.createNode('follicle', name = follicle_name)
        follicleTransform = cmds.listRelatives(follicleShp, parent = True)[0]
        cmds.connectAttr(loft[0] + '.worldMatrix[0]', follicleShp + '.inputWorldMatrix')
        cmds.connectAttr(loft[0] + '.local', follicleShp + '.inputSurface')
        cmds.setAttr(follicleShp + ".parameterU", 0.5)
        cmds.setAttr(follicleShp + ".parameterV", pV_num)
        pV_num += 1.0 / num_edges if num_edges else 0  # Prevent division by zero
        cmds.connectAttr(follicleShp + '.outRotate', follicleTransform + '.rotate')
        cmds.connectAttr(follicleShp + '.outTranslate', follicleTransform + '.translate')
        cmds.setAttr(follicleShp + '.simulationMethod', 2)
        cmds.parent(follicleShp, fol_grp)

        #Create group
        fol_ctrl_grp = cmds.group(em = True, name = fol_ctrl_grp_name)
        cmds.parent(fol_ctrl_grp, follicleShp)

        #Create joint
        joint_name = cmds.joint(radius = 0.2, name = joint_name)
        
        #cmds.parent(joint_name, world = True)
        fol_joint_position = cmds.xform(follicleTransform, query = True, worldSpace = True, translation = True)
        cmds.xform(fol_ctrl_grp, worldSpace = True, translation = fol_joint_position)
        cmds.hide(joint_name)
        joints.append(joint_name)

    cmds.parent(fol_grp, ribbon_noTransform_grp)
 
    jnt_grp = cmds.group(em = True, n = "Jnt_grp")

    root_offset_ctrl_group = cmds.group(em = True, n = "Root_offset_grp")

    all_follicles = cmds.listRelatives(fol_grp, c = True) or []

    root_controller = cmds.circle(r = 5, name = "Root_ctrl")[0]
    cmds.setAttr(root_controller + ".rotateY", 90)
    cmds.makeIdentity( apply = True, rotate = True)

    cmds.parent(root_controller, root_offset_ctrl_group)

    loc_grp = cmds.group(em = True, n = "Loc_grp")
    cmds.parent(loc_grp, ribbon_noTransform_grp)

    for index, follicle in enumerate(all_follicles):
        mid_index = len(all_follicles)//2
        if index == 0 :
            controller_joint_position = cmds.getAttr(follicle + ".translate")[0]
            if controller_joint_position not in created_circle_root:
                cmds.xform(root_controller, worldSpace = True, translation = controller_joint_position)
                cmds.makeIdentity(root_controller, apply = True, translate = True, rotate = True, scale = True)
            createCtrl("start", cmds.getAttr(follicle + ".translate")[0], jnt_grp, root_controller, "start_loc", ribbon_noTransform_grp, index, ctrl_dict, ctrl_constraint_dict)
            cmds.parent("start_loc", loc_grp)
            cmds.hide("start_loc")

        elif index == mid_index :
            createCtrl("mid", cmds.getAttr(follicle + ".translate")[0], jnt_grp, root_controller, "mid_loc", ribbon_noTransform_grp, index, ctrl_dict, ctrl_constraint_dict) 
            cmds.parent("mid_loc", loc_grp)
            cmds.hide("mid_loc")

        elif index == len(all_follicles) - 1 :
            createCtrl("end", cmds.getAttr(follicle + ".translate")[0], jnt_grp, root_controller, "end_loc", ribbon_noTransform_grp, index, ctrl_dict, ctrl_constraint_dict)
            cmds.parent("end_loc", loc_grp)
            cmds.hide("end_loc")

    num_elements_each_side = (len(all_follicles)-1)//2
    mid_upper = ((len(all_follicles) - 1) + mid_index) // 2
    mid_lower = (1 + mid_index)//2 

    for i in range(2, num_elements_each_side -2):
        index_right = mid_index + i
        index_left = mid_index - i

        if index_left == mid_lower:
            createCtrl("A", cmds.getAttr(all_follicles[index_left] + ".translate")[0], jnt_grp, root_controller, "A_loc", ribbon_noTransform_grp, index_left, ctrl_dict, ctrl_constraint_dict)
            cmds.parent("A_loc", loc_grp)
            cmds.hide("A_loc")

        if index_right == mid_upper:
            createCtrl("B", cmds.getAttr(all_follicles[index_right] + ".translate")[0], jnt_grp, root_controller, "B_loc", ribbon_noTransform_grp, index_right, ctrl_dict, ctrl_constraint_dict)
            cmds.parent("B_loc", loc_grp)
            cmds.hide("B_loc")

    cmds.parent(jnt_grp, ribbon_transform_grp)
    cmds.parent(root_offset_ctrl_group, ribbon_transform_grp)

    joint_control_created = cmds.listRelatives(jnt_grp, c = True)

    for index, jnt in enumerate(joint_control_created):
        circle_parent_constraint = cmds.listConnections(jnt , type = "parentConstraint")
        if not circle_parent_constraint and cmds.listRelatives(jnt, parent = True):
            cmds.parent(jnt, world = True)
            cmds.delete(jnt)

    joint_in_group = cmds.ls(cmds.listRelatives(jnt_grp, type = "joint"), type = "joint")
    cmds.select(joint_in_group,loft)
    loft_skinCluster = cmds.skinCluster(n = "loft_skinCluster", tsb = True, bindMethod = 0)[0]

    #joint_fol_group = cmds.ls(cmds.listRelatives(fol_jnt_group, type = "joint"), type = "joint")
    cmds.select(joints, select_geo)
    geo_skin = cmds.skinCluster(n = "geo_skinCluster", tsb = True, bindMethod = 0)[0]

    if addFk == True:
        makeFkCtrl(ctrl_dict, ctrl_constraint_dict,root_controller)

    if addDeformation == True:
        makeDeformation(root_controller, loft, ribbon_noTransform_grp, loc_grp, num_edges, all_follicles, fol_grp, ribbon_transform_grp)
        cmds.reorderDeformers("loft_skinCluster", "ribbon_general_BS", "loft_created")


def createCtrl(name_prefix, target_position, parent_group, root_controller, loc_prefix, ribbon_noTransform_grp, index, ctrl_dict, ctrl_constraint_dict):
    created_circle = {}
    created_joint = {}

    ctrl_offset_grp = cmds.group(em = True, n = name_prefix + "_ctrl_offset_grp")
    ctrl_constraint = cmds.group(em = True, n = name_prefix + "_ctrl_const_grp")
    
    main_controller = cmds.circle(r = 1, name = name_prefix + "_ctrl")[0]
    cmds.setAttr(main_controller + ".rotateY", 90)
    cmds.makeIdentity(apply = True, rotate = True)
    cmds.xform(main_controller, centerPivots = True)
    cmds.delete(main_controller, ch = True)
    cmds.parent(main_controller, ctrl_constraint)
    cmds.parent(ctrl_constraint, ctrl_offset_grp)
    
    main_joint = cmds.joint(name = name_prefix + "_joint")
    cmds.xform(main_joint, centerPivots = True)
    cmds.delete(main_joint, ch = True)
    cmds.parent(main_joint, parent_group)

    if target_position not in created_circle:
        cmds.xform(main_controller, worldSpace = True, translation = target_position)
        cmds.makeIdentity(main_controller, apply = True, translate = True, rotate = True, scale = True)
        cmds.xform(main_joint, worldSpace = True, translation = target_position)
        created_circle[target_position] = main_controller
        created_joint[target_position] = main_joint
        cmds.parentConstraint(main_controller, main_joint, mo = True, w = 1)

    cmds.parent(ctrl_offset_grp, root_controller)
    loc = cmds.spaceLocator(name = loc_prefix)
    parentConst = cmds.parentConstraint(main_controller, loc, weight = 1)
    cmds.delete(parentConst)

    ctrl_dict[index] = main_controller
    ctrl_constraint_dict[index] = ctrl_constraint
    
    return ctrl_dict

def makeFkCtrl(ctrl_dict, ctrl_constraint_dict, root_controller):
    parent_offset_grp = cmds.group(empty=True, name="FK_Ctrl_offset_grp")
    previous_ctrl = None
    duplicated_ctrl_dict = {}

    for index, key in enumerate(sorted(ctrl_dict.keys())):
        sorted_controllers = ctrl_dict[key]
        if not sorted_controllers:
            print("No objects found in sorted_controllers.")
            continue
        
        duplicate_ctrl = cmds.duplicate(sorted_controllers)
        if not duplicate_ctrl:
            print("Failed to duplicate controller:", sorted_controller)
            continue

        cmds.setAttr(duplicate_ctrl[0] + ".scaleX", 2)
        cmds.setAttr(duplicate_ctrl[0] + ".scaleY", 2)
        cmds.setAttr(duplicate_ctrl[0] + ".scaleZ", 2)

        new_name = cmds.rename(duplicate_ctrl[0], "Fk_" + sorted_controllers)
        cmds.parent(new_name, world = True)

        duplicated_ctrl_dict[index] = new_name

    for index, key in enumerate(sorted(duplicated_ctrl_dict.keys())):
        dup_ctrl = duplicated_ctrl_dict[key]

        new_fk_offset_grp = cmds.group(empty=True, name="Fk_" + dup_ctrl + "_offset_grp")
        
        if index  == 0 :
            cmds.parent(dup_ctrl, new_fk_offset_grp)
            cmds.parent(new_fk_offset_grp, parent_offset_grp)
            cmds.parent(parent_offset_grp, root_controller)
        else:
            previous_key = list(duplicated_ctrl_dict.keys())[index - 1]
            previous_ctrl = duplicated_ctrl_dict[previous_key]
            cmds.parent(dup_ctrl, new_fk_offset_grp)
            cmds.parent(new_fk_offset_grp, previous_ctrl)

    for index, key in enumerate(sorted(ctrl_constraint_dict.keys())):
        sorted_controllers = ctrl_constraint_dict[key]

    for i, (ctrl_key, constraint_key) in enumerate(zip(sorted(duplicated_ctrl_dict.keys()), sorted(ctrl_constraint_dict.keys()))):
        ctrl = duplicated_ctrl_dict[ctrl_key]
        constraint = ctrl_constraint_dict[constraint_key]
        cmds.parentConstraint(ctrl, constraint, mo=True, w=1)



def makeDeformation(root_controller, loft, ribbon_noTransform_grp, loc_grp, num_edges, all_follicles, fol_grp, ribbon_transform_grp):
    cmds.addAttr(root_controller, longName = "TwistSep", niceName = "------", attributeType = "enum", enumName = "Twist:")
    
    cmds.setAttr(root_controller + ".TwistSep", keyable=True)
    cmds.setAttr(root_controller + ".TwistSep", lock = True)

    cmds.addAttr(root_controller, longName = "Twist", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".Twist", keyable = True)

    cmds.addAttr(root_controller, longName = "TwistOffset", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".TwistOffset", keyable = True)

    cmds.addAttr(root_controller, longName = "affectMid", attributeType = "double", defaultValue = 10)
    cmds.setAttr(root_controller + ".affectMid", keyable = True)

    cmds.addAttr(root_controller, longName = "endTwist", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".endTwist", keyable = True)

    cmds.addAttr(root_controller, longName = "twistAmplitude", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".twistAmplitude", keyable = True)

    cmds.addAttr(root_controller, longName = "twistPosition", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".twistPosition", keyable = True)

    cmds.addAttr(root_controller, longName = "Roll", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".Roll", keyable = True)

    cmds.addAttr(root_controller, longName = "RollOffset", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".RollOffset", keyable = True)

    #cmds.setAttr(root_controller + ".RollOffset", lock = True, keyable = False, channelBox = False)

    cmds.addAttr(root_controller, longName = "RollSep", niceName = "-------", attributeType = "enum", enumName = "Roll:")
    cmds.setAttr(root_controller + ".RollSep", keyable = True)

    cmds.setAttr(root_controller + ".RollSep", lock = True)

    cmds.addAttr(root_controller, longName = "RollAmplitude", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".RollAmplitude", keyable = True)

    cmds.addAttr(root_controller, longName = "RollPosition", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".RollPosition", keyable = True)

    cmds.addAttr(root_controller, longName = "RollAngle", attributeType = "double", minValue = -10, maxValue = 10, defaultValue = 0)
    cmds.setAttr(root_controller + ".RollAngle", keyable = True)

    cmds.addAttr(root_controller, longName = "RollScale", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".RollScale", keyable = True)
    
    cmds.addAttr(root_controller, longName = "RollHeight", attributeType = "double", defaultValue = 1)
    cmds.setAttr(root_controller + ".RollHeight", keyable = True)

    cmds.addAttr(root_controller, longName = "RollTwist", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".RollTwist", keyable = True)

    cmds.addAttr(root_controller, longName = "RollStartDropoff", attributeType = "double", minValue = -10, maxValue = 0, defaultValue = -1)
    cmds.setAttr(root_controller + ".RollStartDropoff", keyable = True)

    cmds.addAttr(root_controller, longName = "RollEndDropoff", attributeType = "double", minValue = 0, maxValue = 10, defaultValue = 1)
    cmds.setAttr(root_controller + ".RollEndDropoff", keyable = True)

    cmds.addAttr(root_controller, longName = "SineSep", niceName = "-------", attributeType = "enum", enumName = "Sine:")
    cmds.setAttr(root_controller + ".SineSep", keyable = True)

    cmds.setAttr(root_controller + ".SineSep", lock = True)

    cmds.addAttr(root_controller, longName = "Amplitude", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".Amplitude", keyable = True)

    cmds.addAttr(root_controller, longName = "SineOffset", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".SineOffset", keyable = True)

    cmds.addAttr(root_controller, longName = "SineTwist", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".SineTwist", keyable = True)
    
    cmds.addAttr(root_controller, longName = "SinePosition", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".SinePosition", keyable = True)

    cmds.addAttr(root_controller, longName = "SineDropoff", attributeType = "double", minValue = -1, maxValue = 1, defaultValue = 1)
    cmds.setAttr(root_controller + ".SineDropoff", keyable = True)

    cmds.addAttr(root_controller, longName = "SineLength", attributeType = "double", minValue = 0, defaultValue = 2)
    cmds.setAttr(root_controller + ".SineLength", keyable = True)

    cmds.addAttr(root_controller, longName = "SineScale", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".SineScale", keyable = True)

    cmds.addAttr(root_controller, longName = "VolumeSep", niceName = "-------", attributeType = "enum", enumName = "Volume:")
    cmds.setAttr(root_controller + ".VolumeSep", keyable = True)

    cmds.setAttr(root_controller + ".VolumeSep", lock = True)

    cmds.addAttr(root_controller, longName = "Volume", attributeType = "double", minValue = -2, maxValue = 2, defaultValue = 0)
    cmds.setAttr(root_controller + ".Volume", keyable = True)

    cmds.addAttr(root_controller, longName = "VolumeMultipler", attributeType = "double", minValue = 0, maxValue = 10, defaultValue = 0)
    cmds.setAttr(root_controller + ".VolumeMultipler", keyable = True)

    cmds.addAttr(root_controller, longName = "StartDropoff", attributeType = "double", minValue = 0, maxValue = 1, defaultValue = 0)
    cmds.setAttr(root_controller + ".StartDropoff", keyable = True)

    cmds.addAttr(root_controller, longName = "EndDropoff", attributeType = "double", minValue = 0, maxValue = 1, defaultValue = 0)
    cmds.setAttr(root_controller + ".EndDropoff", keyable = True)

    cmds.addAttr(root_controller, longName = "VolumeScale", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".VolumeScale", keyable = True)

    cmds.addAttr(root_controller, longName = "VolumePosition", attributeType = "double", defaultValue = 0)
    cmds.setAttr(root_controller + ".VolumePosition", keyable = True)


    ribbon_twist_geo = cmds.duplicate(loft, rr = True, n = "ribbon_twist_geo")

    deformer_geo_grp = cmds.group(em = True, n = "geo_grp")
    cmds.parent(ribbon_twist_geo, deformer_geo_grp)
    cmds.parent(deformer_geo_grp, ribbon_noTransform_grp)
    cmds.hide(ribbon_twist_geo)

    cmds.select(ribbon_twist_geo, replace = True)
    twist_deformer = cmds.nonLinear(type = "twist")[0]
    cmds.setAttr("twist1Handle.rotateZ", 90)

    deformer_grp = cmds.group(em = True, n = "Deformer_grp")
    twist_offset_grp = cmds.group(empty = True, name = "twist_offset_grp")
    cmds.parent("twist1Handle", twist_offset_grp)
    cmds.parent(twist_offset_grp, deformer_grp)
    cmds.parent(deformer_grp, ribbon_noTransform_grp)
    cmds.hide("twist1Handle")

    ribbon_twist_start_sum_pma = cmds.createNode("plusMinusAverage", name = "ribbon_twist_start_sum_pma")
    cmds.connectAttr(root_controller + ".Twist", ribbon_twist_start_sum_pma + ".input1D[0]")
    cmds.connectAttr(root_controller + ".TwistOffset", ribbon_twist_start_sum_pma + ".input1D[1]")
    cmds.connectAttr(root_controller + ".Roll", ribbon_twist_start_sum_pma + ".input1D[2]")
    cmds.connectAttr(root_controller + ".RollOffset", ribbon_twist_start_sum_pma + ".input1D[3]")

    ribbon_twist_end_sum_pma = cmds.createNode("plusMinusAverage", name = "ribbon_twist_end_sum_pma")

    cmds.connectAttr(root_controller + ".endTwist", ribbon_twist_end_sum_pma +".input1D[0]")

    cmds.connectAttr(ribbon_twist_start_sum_pma + ".output1D", "twist1.startAngle")
    cmds.connectAttr(ribbon_twist_end_sum_pma + ".output1D", "twist1.endAngle")

    
    affectMid_mdl = cmds.shadingNode('multDoubleLinear', asUtility = 1, name = 'twist_top_affect_mdl')
    cmds.setAttr(affectMid_mdl + ".input1", -0.1)
    cmds.connectAttr(root_controller + ".affectMid", affectMid_mdl + ".input2")
    cmds.connectAttr(affectMid_mdl + ".output", "twist1.lowBound")

    
    endAffectMid_mdl = cmds.shadingNode('multDoubleLinear', asUtility = 1, name = "twist_end_affect_mdl")
    cmds.setAttr(endAffectMid_mdl + '.input1', 0.1)
    cmds.connectAttr(root_controller + ".affectMid", endAffectMid_mdl + '.input2')
    cmds.connectAttr(endAffectMid_mdl + '.output', 'twist1.highBound')
    
    
    
    sumAmplitudePma = cmds.shadingNode('plusMinusAverage', asUtility = 1, name = 'twist_amplitude_pma')
    cmds.connectAttr(root_controller + ".twistAmplitude", sumAmplitudePma + '.input1D[0]')
    TwistScale = cmds.getAttr("twist1Handle.sy")
    cmds.setAttr(sumAmplitudePma + '.input1D[1]', TwistScale)
    cmds.connectAttr(sumAmplitudePma + '.output1D', "twist1Handle.sy")
    
    cmds.connectAttr(root_controller + ".twistPosition", twist_offset_grp + ".translateX")

    #------------------------------------------------------ Roll ---------------------------------------------------------

    ribbon_roll_geo = cmds.duplicate(loft, rr = True, n = "ribbon_roll_geo")
    cmds.parent(ribbon_roll_geo, deformer_geo_grp)
    cmds.hide(ribbon_roll_geo)

    cmds.select(ribbon_roll_geo, replace = True)
    twist_deformer = cmds.nonLinear(type = "bend")[0]
    cmds.setAttr("bend1Handle.rotateZ", 90)

    bend_offset_grp = cmds.group(empty = True, name = "bend_offset_grp")
    cmds.parent("bend1Handle", bend_offset_grp)
    cmds.parent(bend_offset_grp, deformer_grp)
    cmds.hide("bend1Handle")

    cmds.connectAttr(root_controller + ".RollAmplitude",  'bend1.curvature')
    cmds.connectAttr(root_controller + ".RollStartDropoff",  'bend1.lowBound')
    cmds.connectAttr(root_controller + ".RollEndDropoff",  'bend1.highBound')

    
    cmds.connectAttr(root_controller + ".RollAngle",  bend_offset_grp + '.rotateZ')
    cmds.connectAttr(root_controller + ".RollTwist",  bend_offset_grp + '.rotateY')
    cmds.connectAttr(root_controller + ".RollScale",  bend_offset_grp + '.translateY')
    cmds.connectAttr(root_controller + ".RollHeight",  bend_offset_grp + '.scaleX')
    cmds.connectAttr(root_controller + ".RollPosition",  bend_offset_grp + '.translateX')
    
    #------------------------------------------------------ Sine ---------------------------------------------------------

    ribbon_sine_geo = cmds.duplicate(loft, rr = True, n = "ribbon_sine_geo")
    cmds.parent(ribbon_sine_geo, deformer_geo_grp)
    cmds.hide(ribbon_sine_geo)

    cmds.select(ribbon_sine_geo, replace = True)
    twist_deformer = cmds.nonLinear(type = "sine")[0]
    cmds.setAttr("sine1Handle.rotateZ", 90)

    sine_offset_grp = cmds.group(empty = True, name = "sine_offset_grp")
    cmds.parent("sine1Handle", sine_offset_grp)
    cmds.parent(sine_offset_grp, deformer_grp)
    cmds.hide("sine1Handle")

    sumSineScalePma = cmds.shadingNode('plusMinusAverage', asUtility = 1, name ='sine_scale_pma')
    cmds.connectAttr(root_controller + ".SineScale", sumSineScalePma + '.input1D[0]')
    TwistScale = cmds.getAttr("sine1Handle.sy")
    cmds.setAttr(sumSineScalePma + '.input1D[1]', TwistScale)
    cmds.connectAttr(sumSineScalePma + '.output1D', "sine1Handle.sy")
    
    cmds.connectAttr(root_controller + ".Amplitude", "sine1.amplitude")
    cmds.connectAttr(root_controller + ".SineOffset", "sine1.offset")
    cmds.connectAttr(root_controller + ".SineTwist", "sine1Handle.rotateY")
    cmds.connectAttr(root_controller + ".SineLength", "sine1.wavelength")
    cmds.connectAttr(root_controller + ".SineDropoff", "sine1.dropoff")
    cmds.connectAttr(root_controller + ".SinePosition", sine_offset_grp +".translateX")


    pos1 = cmds.xform("start_loc", q = True, ws = True, t = True)
    pos2 = cmds.xform("mid_loc", q = True, ws = True, t = True)
    pos3 = cmds.xform("end_loc", q = True, ws = True, t = True)


    #Curve'ün üstünden 2 vertex noktası seçip culuster oluştur.ClusterShape'in origin değerini ve pivot noktasını değiştirme.

    curve = cmds.curve(d = 2, p = [(pos1[0], pos1[1], pos1[2]), (pos2[0], pos2[1], pos2[2]), (pos3[0], pos3[1], pos3[2])], n = "ribbon_wire_crv")
    cmds.hide(curve)

    curve_grp = cmds.group(em = True, n = "ribbon_msc_grp")
    cmds.parent(curve, curve_grp)
    cmds.parent(curve_grp,ribbon_noTransform_grp)

    cluster_grp = cmds.group(em = True, n = "cluster_grp")

    ribbon_end_cluster = cmds.cluster(curve + ".cv[1:2]", n = "ribbon_end_cluster")
    cmds.parent(ribbon_end_cluster, cluster_grp)
    cmds.parent(cluster_grp, ribbon_noTransform_grp)
    cmds.hide(ribbon_end_cluster)

    trX_end_loc = cmds.getAttr("end_loc.translateX")

    cmds.setAttr("ribbon_end_clusterHandleShape.originX", trX_end_loc)

    current_pivot_point = cmds.xform("end_loc", q = True, ws = True, t = True)

    cmds.move(current_pivot_point[0], current_pivot_point[1], current_pivot_point[2], 'ribbon_end_clusterHandle.scalePivot', 'ribbon_end_clusterHandle.rotatePivot', absolute=True)

    cmds.connectAttr("end_ctrl.translate", "ribbon_end_clusterHandle.translate")

    #-------------------------------------

    ribbon_start_cluster = cmds.cluster(curve + ".cv[0:1]", n = "ribbon_start_cluster")
    cmds.parent(ribbon_start_cluster, cluster_grp)
    cmds.hide(ribbon_start_cluster)

    trX_start_loc = cmds.getAttr("start_loc.translateX")

    cmds.setAttr("ribbon_start_clusterHandleShape.originX", trX_start_loc)

    current_pivot_point = cmds.xform("start_loc", q = True, ws = True, t = True)

    cmds.move(current_pivot_point[0], current_pivot_point[1], current_pivot_point[2], 'ribbon_start_clusterHandle.scalePivot', 'ribbon_start_clusterHandle.rotatePivot', absolute=True)

    cmds.connectAttr("start_ctrl.translate", "ribbon_start_clusterHandle.translate")

    #-------------------------------------

    ribbon_mid_cluster = cmds.cluster(curve +".cv[1]", n = "ribbon_mid_cluster")
    cmds.parent(ribbon_mid_cluster, cluster_grp)
    cmds.hide(ribbon_mid_cluster)

    cmds.connectAttr("mid_ctrl.translate", "ribbon_mid_clusterHandle.translate")

    cmds.percent('ribbon_end_cluster', 'ribbon_wire_crv.cv[1]', v = 0.5)

    ribbon_wire_geo = cmds.duplicate(loft, rr = True, n = "ribbon_wire_geo")
    cmds.parent(ribbon_wire_geo, deformer_geo_grp)
    cmds.hide(ribbon_wire_geo)
 
    cmds.select(ribbon_wire_geo, curve)

    cmds.wire(ribbon_wire_geo, gw = False, en = 1, ce = 0, li = 0, w = curve)
    cmds.setAttr("wire1.dropoffDistance[0]", 10)

    #------------------------------------------------------ Volume ---------------------------------------------------------

    ribbon_squash_geo = cmds.duplicate(loft, rr = True, n = "ribbon_squash_geo")
    cmds.parent(ribbon_squash_geo, deformer_geo_grp)
    cmds.hide(ribbon_squash_geo)

    fol_squash_grp = cmds.group(em = True, name = "follicle_squash_grp")
    pV_squash_num = 0.0
    cmds.parent(fol_squash_grp, ribbon_noTransform_grp)
    transZ_value = 0

    ctrl_fol_grp = []

    for follicle_ctrl in all_follicles:
        follicle_ctrl_grp = cmds.listRelatives(follicle_ctrl, c=True) or []
        for x in follicle_ctrl_grp:
            print(x)

        ctrl_fol_grp.append(x)

    for num in range(num_edges + 1):

        # Construct unique names
        follicle_squash_name = "follicle_squash_{0}".format(num)
        
        follicleSquashShp = cmds.createNode('follicle', name = follicle_squash_name)
        follicleSquashTransform = cmds.listRelatives(follicleSquashShp, parent = True)[0]
        cmds.connectAttr(ribbon_squash_geo[0] + '.worldMatrix[0]', follicleSquashShp + '.inputWorldMatrix')
        cmds.connectAttr(ribbon_squash_geo[0] + '.local', follicleSquashShp + '.inputSurface')
        cmds.setAttr(follicleSquashShp + ".parameterU", 0)
        cmds.setAttr(follicleSquashShp + ".parameterV", pV_squash_num)
        pV_squash_num += 1.0 / num_edges if num_edges else 0  # Prevent division by zero
        cmds.connectAttr(follicleSquashShp + '.outRotate', follicleSquashTransform + '.rotate')
        cmds.connectAttr(follicleSquashShp + '.outTranslate', follicleSquashTransform + '.translate')
        cmds.setAttr(follicleSquashShp + '.simulationMethod', 2)
        cmds.parent(follicleSquashShp, fol_squash_grp)

    cmds.hide(fol_squash_grp)

    cmds.select(ribbon_squash_geo, replace = True)
    squash_deformer = cmds.nonLinear(type = "squash")[0]
    cmds.setAttr("squash1Handle.rotateZ", 90)
    squash_offset_grp = cmds.group(empty = True, name = "squash_offset_grp")
    cmds.parent("squash1Handle", squash_offset_grp)
    cmds.parent(squash_offset_grp, deformer_grp)
    cmds.hide("squash1Handle")

    #cmds.connectAttr("Mid_Ctrl.Volume", 'squash1.factor')
    cmds.connectAttr(root_controller + ".StartDropoff", 'squash1.startSmoothness')
    cmds.connectAttr(root_controller + ".EndDropoff", 'squash1.endSmoothness')

    cmds.connectAttr(root_controller + ".VolumePosition", squash_offset_grp + ".translateX")

    all_follicles_for_connection = cmds.listRelatives(fol_squash_grp, c=True) or []

    for index, follicle in enumerate(all_follicles_for_connection):
        mult_squash = cmds.shadingNode("multDoubleLinear", asUtility = True, name = follicle + "_squashMultiplier_md")
        pma_squash = cmds.shadingNode("plusMinusAverage", asUtility = True, name = follicle + "_squashSum_pma")

        #mult_squash_mtd = cmds.shadingNode("multiplyDivide", asUtility=True, name = follicle + "_squashMultiplier_mtd")
        cmds.connectAttr(root_controller + ".VolumeMultipler", mult_squash + ".input1")
        cmds.connectAttr(follicle + ".translateZ", mult_squash + ".input2")

        #cmds.setAttr(mult_squash_mtd + ".operation", 2)
        cmds.setAttr(pma_squash + ".input1D[0]", 1)
        cmds.connectAttr(mult_squash + ".output", pma_squash + ".input1D[1]")

        index = 0

        while index < len(ctrl_fol_grp):
            ctrl_grp_index = ctrl_fol_grp[index]
            # Check if scaleY and scaleZ attributes have incoming connections
            if not cmds.listConnections(ctrl_grp_index + ".scaleY", source=True) \
                and not cmds.listConnections(ctrl_grp_index + ".scaleZ", source=True):
            # Connect to scaleY and scaleZ attributes if they have no incoming connections
                    cmds.connectAttr(pma_squash + ".output1D", ctrl_grp_index + ".scaleY")
                    cmds.connectAttr(pma_squash + ".output1D", ctrl_grp_index + ".scaleZ")
                    
                    break
            index += 1

    squash_reverse_mdl = cmds.shadingNode("multDoubleLinear", asUtility=True, name = "ribbon_squash_reverse_mdl")
    cmds.connectAttr(root_controller + ".Volume", squash_reverse_mdl + ".input1" )
    cmds.setAttr(squash_reverse_mdl + ".input2", -1)

    cmds.connectAttr(squash_reverse_mdl + ".output", "squash1.factor")

    pma_scale = cmds.shadingNode("plusMinusAverage", asUtility=True, name = "ribbon_scale_squash_sum_pma")
    cmds.connectAttr(root_controller + ".VolumeScale", "ribbon_scale_squash_sum_pma.input1D[0]")
    scale_value = cmds.getAttr("squash1Handle.scaleX")
    print(scale_value)
    cmds.setAttr("ribbon_scale_squash_sum_pma.input1D[1]", scale_value)
    cmds.connectAttr("ribbon_scale_squash_sum_pma.output1D", "squash1Handle.scaleY")

    all_follicles_for_scale = cmds.listRelatives(fol_grp, c = True) or []

    for index, follicle in enumerate(all_follicles_for_scale):

        cmds.scaleConstraint(ribbon_transform_grp, all_follicles_for_scale[index], offset = (1,1,1), weight = 1,)

    cmds.blendShape(ribbon_twist_geo, ribbon_roll_geo, ribbon_sine_geo, ribbon_wire_geo, ribbon_squash_geo, loft, n = "ribbon_general_BS", en=1, tc=True, w=[(0, 1.0)])
    cmds.setAttr("ribbon_general_BS.ribbon_twist_geo", 1)
    cmds.setAttr("ribbon_general_BS.ribbon_roll_geo", 1)
    cmds.setAttr("ribbon_general_BS.ribbon_sine_geo", 1)
    cmds.setAttr("ribbon_general_BS.ribbon_wire_geo", 1)
    cmds.setAttr("ribbon_general_BS.ribbon_squash_geo", 1)

def process_textField(text_field, text_field_02, curve_group_name, loft_name):
    selected_edges = cmds.textField(text_field, query = True, text = True)
    select_geo = cmds.textField(text_field_02, query = True, text = True)
    
    selected_edges_list = selected_edges.split(',')
    cmds.select(selected_edges_list)

    if not selected_edges:
        cmds.warning(f"No edges selected in {text_field}.")
        return None, None

    num_edges = len(selected_edges.split(','))
    selected_edges = cmds.ls(selection = True)

    curve = cmds.polyToCurve(form=2, degree=3)[0]
    cmds.xform(curve, centerPivots = True)
    cmds.delete(curve, ch = True)

    curve_grp = cmds.group(em = True, name = curve_group_name)
    cmds.parent(curve, curve_grp)
    cmds.hide(curve)

    dup_curve = cmds.duplicate(curve, smartTransform = True)[0]
    cmds.xform(dup_curve, centerPivots = True)
    cmds.delete(dup_curve, ch = True)

    cmds.setAttr(dup_curve + ".translateZ", 0.1)
    loft = cmds.loft(curve, dup_curve, n = loft_name, ch = True, polygon = 0)[0]
    createFol(num_edges, [loft],select_geo)

    return selected_edges, num_edges

def addFk(state):
    global addFk
    addFk = state

def addDeformation(state):
    global addDeformation
    addDeformation = state

def build(*args):
    first_selected, first_num_edges = process_textField(first_tool, model_tool, "curve_grp", "loft_created")
    if not first_selected:
        cmds.warning("No edges selected in either text field.") 

def window_setup():
    global model_tool, first_tool, addFk_checkbox, addDeformation

    if cmds.window("RibbonTool", exists = True):
        cmds.deleteUI("RibbonTool")

    window = cmds.window("RibbonTool", title = "RibbonTool", widthHeight = (300,100))
    cmds.columnLayout(adjustableColumn = True)

    cmds.rowColumnLayout(numberOfColumns = 3, columnWidth = [(1, 100), (2, 150), (3, 100)], columnSpacing = [(1, 10), (2, 10), (3, 10)], rowSpacing = [1, 5])
    cmds.text(label = "Model:")
    model_tool = cmds.textField()
    cmds.button(label = "Select Model", command = select_model)
    cmds.setParent("..")

    cmds.rowColumnLayout(numberOfColumns = 3, columnWidth = [(1, 100), (2, 150), (3, 100)], columnSpacing = [(1, 10), (2, 10), (3, 10)], rowSpacing = [1, 5])
    cmds.text(label = "Edges:")
    first_tool = cmds.textField()
    cmds.button(label = "Select Edges", command = select_edges)
    cmds.setParent("..")

    cmds.rowColumnLayout(numberOfColumns = 3, columnWidth = [(1, 100), (2, 150), (3, 150)], columnSpacing = [(1, 10), (2, 10), (3, 10)], rowSpacing = [1, 10])
    cmds.text(label = "")
    addFk_checkbox = cmds.checkBox(label="Add FK", changeCommand = addFk)
    addDeformation = cmds.checkBox(label = "Add Deformation", changeCommand = addDeformation)
    cmds.setParent("..")

    cmds.rowColumnLayout(numberOfColumns = 3, columnWidth = [(1, 100), (2, 150), (3, 50)], columnSpacing=[(1, 10), (2, 10), (3, 10)], rowSpacing=[10, 10])
    cmds.text(label = "")
    cmds.button(label = "Build", command = build, width = 50)  # Smaller width for the button
    cmds.text(label = "")
    cmds.setParent("..")
    cmds.showWindow("RibbonTool")

window_setup()