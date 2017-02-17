#!/bin/bash

#import productionTAG, skimProductionTAG and datasetFile
source configs/2015-v1.conf

#all input files to process
#DYMadIncl
#mcIdentifier='DYJets_madgraph'
#inputFiles=('root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_10_2_KFJ.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_11_2_PR1.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_12_2_DuM.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_13_2_yHC.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_14_2_bWL.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_15_2_rK5.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_16_2_yJp.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_17_2_Gff.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_18_2_KEW.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_19_2_3FO.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_1_2_7zA.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_20_2_hc0.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_21_2_nur.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_22_1_SSm.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_2_2_89c.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_3_2_yNe.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_4_2_aL0.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_5_2_RDU.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_6_2_wPy.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_7_2_YMG.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_8_2_xNp.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_madgraph_SHv12/DYJets_madgraph_9_2_EMH.root')

#DYAMC
#mcIdentifier='DYJets_amctnlo'
#inputFiles=('root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_100_1_Q2b.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_101_1_cT4.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_102_1_omu.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_103_1_10K.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_104_1_ZmN.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_105_1_oBB.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_106_1_JbU.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_107_1_ESz.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_108_1_fCX.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_109_1_LDB.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_10_1_qQL.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_110_1_MMu.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_111_1_3Kv.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_112_1_pLw.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_113_1_WYj.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_114_1_UOz.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_115_1_PMs.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_116_1_IpD.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_117_1_PFE.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_118_1_nlE.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_119_1_6vw.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_11_1_CaR.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_120_1_fZV.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_121_1_Gt8.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_122_1_Crv.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_123_1_k4Q.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_124_1_q4S.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_125_1_NJt.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_126_1_WXz.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_127_1_g3J.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_128_1_nzg.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_129_1_ioj.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_12_1_Web.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_130_1_O9U.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_131_1_BbM.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_132_1_fNI.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_133_1_2AE.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_134_1_hf4.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_135_1_YDX.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_13_1_v8f.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_14_1_yWy.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_15_1_YCE.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_16_1_3BZ.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_17_1_GWn.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_18_1_MaW.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_19_1_1Uf.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_1_1_vjJ.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_20_1_eZG.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_21_1_wEP.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_22_1_NUa.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_23_1_ICc.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_24_1_nf0.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_25_1_CEe.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_26_1_JJg.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_27_1_YLj.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_28_1_zDZ.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_29_1_f5Y.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_2_1_Prv.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_30_1_OuK.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_31_1_y4L.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_32_1_zAc.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_33_1_E5g.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_34_1_NIV.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_35_1_PTY.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_36_1_lAT.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_37_1_LgT.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_38_1_v4S.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_39_1_yKN.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_3_1_gzn.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_40_1_ZQp.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_41_1_xou.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_42_1_hhU.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_43_1_COh.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_44_1_byh.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_45_1_XR8.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_46_1_JL8.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_47_1_fRl.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_48_1_8AY.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_49_1_aI6.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_4_1_Alq.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_50_1_IPN.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_51_1_oFr.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_52_1_CgW.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_53_1_gta.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_54_1_l4p.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_55_1_DtW.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_56_1_amB.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_57_1_5Zj.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_58_1_XHS.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_59_1_Gm5.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_5_1_JE6.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_60_1_Y07.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_61_1_XPe.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_62_1_Ku5.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_63_1_Er0.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_64_1_TOK.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_65_1_Vpk.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_66_1_6St.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_67_1_CFd.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_68_1_PyY.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_69_1_xSe.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_6_1_HkD.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_70_1_15y.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_71_1_b6N.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_72_1_aeu.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_73_1_ZsR.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_74_1_9QK.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_75_1_Nnz.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_76_1_f2K.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_77_1_8C9.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_78_1_fEH.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_79_1_60W.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_7_1_ocU.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_80_1_ZcN.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_81_1_hQh.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_82_1_7iA.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_83_1_t81.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_84_1_Wq7.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_85_1_BeC.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_86_1_NCc.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_87_1_n5S.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_88_1_Twk.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_89_1_9V9.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_8_1_Nra.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_90_1_XEF.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_91_1_6iu.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_92_1_FyW.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_93_1_qrh.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_94_1_9iN.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_95_1_5pZ.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_96_1_Zob.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_97_1_LZo.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_98_1_m3L.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_99_1_0sS.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/DYJets_amctnlo_SHv12/DYJets_amctnlo_9_1_Q2a.root' )

#DYMadHT100to200
#mcIdentifier=('DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200' 'DYJets_madgraph_ht100to200')
#inputFiles=('root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_10_2_kP2.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_11_1_pDt.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_12_1_XBz.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_13_3_911.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_14_2_T4l.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_15_1_pVb.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_16_1_l9I.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_1_2_TRt.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_2_1_KOp.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_3_1_3BN.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_4_2_z5k.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_5_1_3G8.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_6_1_95A.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_7_1_WNZ.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_8_1_Vh3.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht100to200_SHv12/DYJets_madgraph_ht100to200_9_1_aQu.root')

##DYMadHT200to400
#mcIdentifier=('DYJets_madgraph_ht200to400' 'DYJets_madgraph_ht200to400' 'DYJets_madgraph_ht200to400' 'DYJets_madgraph_ht200to400')
#inputFiles=('root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht200to400_SHv12/DYJets_madgraph_ht200to400_1_1_TCb.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht200to400_SHv12/DYJets_madgraph_ht200to400_2_1_2qZ.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht200to400_SHv12/DYJets_madgraph_ht200to400_3_1_Dby.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht200to400_SHv12/DYJets_madgraph_ht200to400_4_1_UVm.root')

##DYMadHT400to600
#mcIdentifier=('DYJets_madgraph_ht400to600' 'DYJets_madgraph_ht400to600' 'DYJets_madgraph_ht400to600' 'DYJets_madgraph_ht400to600')
#inputFiles=('root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht400to600_SHv12/DYJets_madgraph_ht400to600_1_1_NYr.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht400to600_SHv12/DYJets_madgraph_ht400to600_2_1_263.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht400to600_SHv12/DYJets_madgraph_ht400to600_3_1_vfG.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht400to600_SHv12/DYJets_madgraph_ht400to600_4_1_Mtb.root')

##DYMadHT600toInf
#mcIdentifier=('DYJets_madgraph_ht600toInf' 'DYJets_madgraph_ht600toInf' 'DYJets_madgraph_ht600toInf' 'DYJets_madgraph_ht600toInf')
#inputFiles=('root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht600toInf_SHv12/DYJets_madgraph_ht600toInf_1_1_d7e.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht600toInf_SHv12/DYJets_madgraph_ht600toInf_2_1_2RV.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht600toInf_SHv12/DYJets_madgraph_ht600toInf_3_1_5Uz.root' 'root://eoscms.cern.ch//eos/cms/store/caf/user/shervin/skims/DYJets_madgraph_ht600toInf_SHv12/DYJets_madgraph_ht600toInf_4_1_oeT.root')


#all WREEJJ, only one file per skim
#mcIdentifier=('WRtoEEJJ_800_400' 'WRtoEEJJ_1000_500' 'WRtoEEJJ_1200_600' 'WRtoEEJJ_1400_700' 'WRtoEEJJ_1600_800' 'WRtoEEJJ_1800_1400' 'WRtoEEJJ_2000_1000' 'WRtoEEJJ_2200_1100' 'WRtoEEJJ_2400_1200' 'WRtoEEJJ_2600_1300' 'WRtoEEJJ_2800_1400' 'WRtoEEJJ_3000_1500' 'WRtoEEJJ_3200_1600' 'WRtoEEJJ_3600_1800' 'WRtoEEJJ_3800_1900' 'WRtoEEJJ_4000_2000')
#inputFiles=('root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_800_400_SHv12/WRtoEEJJ_800_400_1_1_EnO.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_1000_500_SHv12/WRtoEEJJ_1000_500_1_1_E1f.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_1200_600_SHv12/WRtoEEJJ_1200_600_1_1_zZN.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_1400_700_SHv12/WRtoEEJJ_1400_700_1_1_PwB.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_1600_800_SHv12/WRtoEEJJ_1600_800_1_1_J7N.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_1800_1400_SHv12/WRtoEEJJ_1800_1400_1_1_561.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_2000_1000_SHv12/WRtoEEJJ_2000_1000_1_1_hmG.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_2200_1100_SHv12/WRtoEEJJ_2200_1100_1_1_wPE.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_2400_1200_SHv12/WRtoEEJJ_2400_1200_1_1_XeT.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_2600_1300_SHv12/WRtoEEJJ_2600_1300_1_1_5Kb.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_2800_1400_SHv12/WRtoEEJJ_2800_1400_1_2_Jo9.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_3000_1500_SHv12/WRtoEEJJ_3000_1500_1_2_2El.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_3200_1600_SHv12/WRtoEEJJ_3200_1600_1_2_ji7.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_3600_1800_SHv12/WRtoEEJJ_3600_1800_1_2_XeQ.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_3800_1900_SHv12/WRtoEEJJ_3800_1900_1_2_0fx.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_4000_2000_SHv12/WRtoEEJJ_4000_2000_1_2_Mjq.root')

#all WRMuMuJJ, only one file per skim
#mcIdentifier=('WRtoMuMuJJ_800_400' 'WRtoMuMuJJ_1000_500' 'WRtoMuMuJJ_1200_600' 'WRtoMuMuJJ_1400_700' 'WRtoMuMuJJ_1600_800' 'WRtoMuMuJJ_1800_1400' 'WRtoMuMuJJ_2000_1000' 'WRtoMuMuJJ_2200_1100' 'WRtoMuMuJJ_2400_1200' 'WRtoMuMuJJ_2600_1300' 'WRtoMuMuJJ_2800_1400' 'WRtoMuMuJJ_3000_1500' 'WRtoMuMuJJ_3200_1600' 'WRtoMuMuJJ_3600_1800' 'WRtoMuMuJJ_3800_1900' 'WRtoMuMuJJ_4000_2000')
#inputFiles=('root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_800_400_SHv12/WRtoEEJJ_800_400_1_1_EnO.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_1000_500_SHv12/WRtoEEJJ_1000_500_1_1_E1f.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_1200_600_SHv12/WRtoEEJJ_1200_600_1_1_zZN.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_1400_700_SHv12/WRtoEEJJ_1400_700_1_1_PwB.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_1600_800_SHv12/WRtoEEJJ_1600_800_1_1_J7N.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_1800_1400_SHv12/WRtoEEJJ_1800_1400_1_1_561.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_2000_1000_SHv12/WRtoEEJJ_2000_1000_1_1_hmG.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_2200_1100_SHv12/WRtoEEJJ_2200_1100_1_1_wPE.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_2400_1200_SHv12/WRtoEEJJ_2400_1200_1_1_XeT.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_2600_1300_SHv12/WRtoEEJJ_2600_1300_1_1_5Kb.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_2800_1400_SHv12/WRtoEEJJ_2800_1400_1_2_Jo9.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_3000_1500_SHv12/WRtoEEJJ_3000_1500_1_2_2El.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_3200_1600_SHv12/WRtoEEJJ_3200_1600_1_2_ji7.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_3600_1800_SHv12/WRtoEEJJ_3600_1800_1_2_XeQ.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_3800_1900_SHv12/WRtoEEJJ_3800_1900_1_2_0fx.root' 'root://eoscms.cern.ch//eos/cms/store/user/shervin/skims/WRtoEEJJ_4000_2000_SHv12/WRtoEEJJ_4000_2000_1_2_Mjq.root')

##for all datasets in configs/missing_datasets.dat
##first get all the dataset short names, and save them to a vector
mcIdentifier=(`cat $datasetFile | grep -v '#' | awk '{print $1}'`)
eosSkimPath='/eos/cms/store/group/phys_exotica/leptonsPlusJets/WR/skims/'
eosReadingTag='root://eoscms.cern.ch/'

#now loop over all elements in mcIdentifier
for j in ${!mcIdentifier[*]}
do
	#number used to distinguish different jobs processing the same dataset
	startingCount=1

	#now get all the skim files from the dataset linked to mcIdentifier
	inputFiles=(`eos ls $eosSkimPath${mcIdentifier[$j]}$skimProductionTAG/`)

	for i in ${!inputFiles[*]}
	do
		#echo $eosSkimPath${mcIdentifier[$j]}$skimProductionTAG/${inputFiles[$i]}
		#replace NNN by a number and INPTFILE with a file path name from inputFiles in temp_runAnalysis_cfg.py
		eval "cd test"
		eval "sed 's@NNN@$startingCount@g' temp_runAnalysis_cfg.py > tempOne.py"
		eval "sed 's@TAGNAME@${mcIdentifier[$j]}@g' tempOne.py > tempTwo.py"
		eval "sed 's@INPTFILE@$eosReadingTag$eosSkimPath${mcIdentifier[$j]}$skimProductionTAG/${inputFiles[$i]}@g' tempTwo.py > runAnalysis_cfg_${mcIdentifier[$j]}_${startingCount}.py"
		#eval "sed 's@TAGNAME@$mcIdentifier@g' tempOne.py > tempTwo.py"
		#eval "sed 's@INPTFILE@${inputFiles[$i]}@g' tempTwo.py > runAnalysis_cfg_${mcIdentifier}_${startingCount}.py"
		rm tempOne.py tempTwo.py
		eval "cd .."

		#replace NNN by a number in runAnalysisDYMC.sh
		eval "sed 's@NNN@$startingCount@g' runAnalysisDYMC.sh > tempOne.sh"
		eval "sed 's@TAGNAME@${mcIdentifier[$j]}@g' tempOne.sh > runAnalysisDYMC_${mcIdentifier[$j]}_${startingCount}.sh"
		eval "chmod u+x runAnalysisDYMC_${mcIdentifier[$j]}_${startingCount}.sh"
		#eval "sed 's@TAGNAME@$mcIdentifier@g' tempOne.sh > runAnalysisDYMC_${mcIdentifier}_${startingCount}.sh"
		#eval "chmod u+x runAnalysisDYMC_${mcIdentifier}_${startingCount}.sh"
		rm tempOne.sh

		#submit the job to 1nh queue and request at least 2 GB of disk
		#echo "bsub -R 'pool>2000' -q 8nh -J runAnalysisDYMC_${mcIdentifier[$j]}_Part_${startingCount} < runAnalysisDYMC_${mcIdentifier[$j]}_${startingCount}.sh"
		eval "bsub -R 'pool>2000' -q 8nh -J runAnalysisDYMC_${mcIdentifier[$j]}_Part_${startingCount} < runAnalysisDYMC_${mcIdentifier[$j]}_${startingCount}.sh"
		#eval "bsub -R 'pool>2000' -q 8nh -J runAnalysisDYMC_${mcIdentifier}_Part_${startingCount} < runAnalysisDYMC_${mcIdentifier}_${startingCount}.sh"

		let startingCount=startingCount+1
	done

done
