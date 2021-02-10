//
//  Worm.cpp
//  one
//
//  Created by Eduardo Izquierdo on 9/25/15.
//  Copyright © 2015 Eduardo Izquierdo. All rights reserved.
//

#include "Worm.h"

#define HEADSR
#define VNCSR

int nn(int neuronNumber, int unitNumber)
{
    return neuronNumber+((unitNumber-1)*N_neuronsperunit);
}

// The constructor
Worm::Worm(TVector<double> &v,double output)
{
    // Muscles
    m.SetMuscleParams(N_muscles, T_muscle);

    // Nervous system // Ventral cord
    n.SetCircuitSize(N_units*N_neuronsperunit, 4, 4);

    int db, dd, vba, vda, vbp, vdp;
    int ddNext, dbNext, vdaNext, vbaNext;

    for (int u = 1; u <= N_units; u++){
        db = nn(DB,u);
        dd = nn(DD,u);
        vba = nn(VBA,u);
        vbp = nn(VBP,u);
        vda = nn(VDA,u);
        vdp = nn(VDP,u);

        ddNext = nn(DD,u+1);
        dbNext = nn(DB,u+1);
        vdaNext = nn(VDA,u+1);
        vbaNext = nn(VBA,u+1);

        // Bias
        //  B-
        n.SetNeuronBias(db, v(1));
        n.SetNeuronBias(vba, v(1));
        n.SetNeuronBias(vbp, v(1));
        //  D-
        n.SetNeuronBias(dd, v(2));
        n.SetNeuronBias(vda, v(2));
        n.SetNeuronBias(vdp, v(2));

        // Time Constant
        //  B-
        n.SetNeuronTimeConstant(db, v(3));
        n.SetNeuronTimeConstant(vba, v(3));
        n.SetNeuronTimeConstant(vbp, v(3));
        //  D-
        n.SetNeuronTimeConstant(dd, v(4));
        n.SetNeuronTimeConstant(vda, v(4));
        n.SetNeuronTimeConstant(vdp, v(4));

        // Self connections
        //  B-
        n.SetChemicalSynapseWeight(db, db, v(5));
        n.SetChemicalSynapseWeight(vba, vba, v(5));
        n.SetChemicalSynapseWeight(vbp, vbp, v(5));
        //  D-
        n.SetChemicalSynapseWeight(dd, dd, v(6));
        n.SetChemicalSynapseWeight(vda, vda, v(6));
        n.SetChemicalSynapseWeight(vdp, vdp, v(6));

        // Chemical synapses (COLORS AS IN FIGURE)
        // xB -> xD
        n.SetChemicalSynapseWeight(db, dd, v(7));       // Lighter Green
        n.SetChemicalSynapseWeight(vba, vda, v(7));
        n.SetChemicalSynapseWeight(vbp, vdp, v(7));
        // xB -> yD
        n.SetChemicalSynapseWeight(db, vda, v(8));      // Darker Green
        n.SetChemicalSynapseWeight(db, vdp, v(8));
        n.SetChemicalSynapseWeight(vba, dd, v(8)/2);
        n.SetChemicalSynapseWeight(vbp, dd, v(8)/2);
        // xD- -> yD
        n.SetChemicalSynapseWeight(dd, vda, v(9));     // Darker Blue

        // Gap junctions within the unit
        n.SetElectricalSynapseWeight(dd, vda, v(10));   // Light Gray
        n.SetElectricalSynapseWeight(dd, vdp, v(10));

        // D-
        n.SetElectricalSynapseWeight(vda, vdp, v(11));  // Blue
        // B-
        n.SetElectricalSynapseWeight(vba, vbp, v(12));  // Green

        // Gap junctions across units
        if (u < N_units){
            //  D-
            n.SetElectricalSynapseWeight(dd, ddNext, v(11));    // Blue
            n.SetElectricalSynapseWeight(vdp, vdaNext, v(11));
            //  B-
            n.SetElectricalSynapseWeight(db, dbNext, v(12));  // Green
            n.SetElectricalSynapseWeight(vbp, vbaNext, v(12));
            //
            n.SetElectricalSynapseWeight(vbp, dbNext, v(13)); // Darker Gray
        }
    }

    // Stretch receptor
    sr.SetStretchReceptorParams(N_segments, N_stretchrec, v(14), v(28));

    // NMJ Weight
    NMJ_DB = v(15);
    NMJ_VBa = v(15);
    NMJ_VBp = v(15);
    NMJ_DD = v(16);
    NMJ_VDa = v(16);
    NMJ_VDp = v(16);

    // Head Circuit
    h.SetCircuitSize(4, 3, 2);

    // Bias
    h.SetNeuronBias(SMDD, v(17));
    h.SetNeuronBias(SMDV, v(17));
    h.SetNeuronBias(RMDD, v(18));
    h.SetNeuronBias(RMDV, v(18));

    // Time-Constant
    h.SetNeuronTimeConstant(SMDD, v(19));
    h.SetNeuronTimeConstant(SMDV, v(19));
    h.SetNeuronTimeConstant(RMDD, v(20));
    h.SetNeuronTimeConstant(RMDV, v(20));

    // Self-Connection
    h.SetChemicalSynapseWeight(SMDD, SMDD, v(21));
    h.SetChemicalSynapseWeight(SMDV, SMDV, v(21));
    h.SetChemicalSynapseWeight(RMDD, RMDD, v(22));
    h.SetChemicalSynapseWeight(RMDV, RMDV, v(22));

    // Chemical Synxapses (ALL)
    h.SetChemicalSynapseWeight(SMDD, SMDV, v(23));
    h.SetChemicalSynapseWeight(SMDV, SMDD, v(23));

    h.SetChemicalSynapseWeight(SMDD, RMDV, v(24));
    h.SetChemicalSynapseWeight(SMDV, RMDD, v(24));

    h.SetChemicalSynapseWeight(RMDD, RMDV, v(25));
    h.SetChemicalSynapseWeight(RMDV, RMDD, v(25));

    // Gap Junctions
    h.SetElectricalSynapseWeight(SMDD, RMDD, v(26));
    h.SetElectricalSynapseWeight(SMDV, RMDV, v(26));
    h.SetElectricalSynapseWeight(RMDV, RMDD, v(27));

    // NMJ Weights
    NMJ_SMDD = v(29);
    NMJ_SMDV = v(29);
    NMJ_RMDD = v(30);
    NMJ_RMDV = v(30);

    // NMJ Gain
    NMJ_Gain_Map = 0.5;
    NMJ_Gain.SetBounds(1, N_muscles);
    for (int i=1; i<=N_muscles; i++)
    {
        NMJ_Gain(i) = 0.7*(1.0 - (((i-1)*NMJ_Gain_Map)/N_muscles));
    }
}

void Worm::InitializeState(RandomState &rs)
{
    t = 0.0;
    n.RandomizeCircuitState(-0.5, 0.5, rs);
    h.RandomizeCircuitState(-0.5, 0.5, rs);
    b.InitializeBodyState();
    m.InitializeMuscleState();
}

void Worm::HeadStep(double StepSize, double output)
{
    // Update Nervous System
    h.EulerStep(StepSize);

    // Time
    t += StepSize;
}

void Worm::Step(double StepSize, double output)
{
    int mi;
    int mt;
    double ds, vs;

    double dorsalHeadInput, ventralHeadInput, ventralHeadInputA, ventralHeadInputP;

    // Update Body
    b.StepBody(StepSize);

    // Set input to Stretch Receptors from Body
    for(int i = 1; i <= N_segments; ++i){
        ds = (b.DorsalSegmentLength(i) - b.RestingLength(i))/b.RestingLength(i);
        vs =  (b.VentralSegmentLength(i) - b.RestingLength(i))/b.RestingLength(i);
        sr.SetDorsalInput(i, ds);
        sr.SetVentralInput(i, vs);
    }

    // Update Stretch Receptors
    sr.Update();

    // Set input to Nervous System (Head) from Stretch Receptors
#ifdef HEADSR
    if (output == 1){
        h.SetNeuronExternalInput(SMDD, sr.HeadDorsalOutput());    // Average of first
        h.SetNeuronExternalInput(SMDV, sr.HeadVentralOutput());   // to segments
    }
#endif

    // Set input to Nervous System (Ventral Cord) from Stretch Receptors
#ifdef VNCSR
    for (int i = 1; i <= N_units; i++){
        n.SetNeuronExternalInput(nn(DB,i), sr.VCDorsalOutput(i));
        n.SetNeuronExternalInput(nn(VBA,i), sr.VCVentralAOutput(i));
        n.SetNeuronExternalInput(nn(VBP,i), sr.VCVentralPOutput(i));
    }
#endif

    // Update Nervous System
    h.EulerStep(StepSize);
    n.EulerStep(StepSize);

    // Set input to Muscles
    //  Input from the head circuit
    dorsalHeadInput = NMJ_SMDD*h.NeuronOutput(SMDD) + NMJ_RMDV*h.NeuronOutput(RMDD);
    ventralHeadInput = NMJ_SMDV*h.NeuronOutput(SMDV) + NMJ_RMDD*h.NeuronOutput(RMDV);

    for (int i = 1; i <= HeadMotorNeuronMuscles; i++){
        m.SetDorsalMuscleInput(i, NMJ_Gain(i)*dorsalHeadInput);
        m.SetVentralMuscleInput(i, NMJ_Gain(i)*ventralHeadInput);
    }

    // Set input to Muscles from Ventral Cord
    //  Dorsal muscles (each motor neuron innervates three muscles, no overlap)
    for (int i = VNCMuscleStart; i <= N_muscles; i++){
        mi = (int) ((i-VNCMuscleStart)/NmusclePerNU)+1;
        dorsalHeadInput = NMJ_DD*n.NeuronOutput(nn(DD,mi)) + NMJ_DB*n.NeuronOutput(nn(DB,mi));
        m.SetDorsalMuscleInput(i, NMJ_Gain(i)*dorsalHeadInput);
        ventralHeadInputA = NMJ_VDa*n.NeuronOutput(nn(VDA,mi)) + NMJ_VBa*n.NeuronOutput(nn(VBA,mi));
        ventralHeadInputP = NMJ_VDp*n.NeuronOutput(nn(VDP,mi)) + NMJ_VBp*n.NeuronOutput(nn(VBP,mi));
        mt = (i-VNCMuscleStart)%NmusclePerNU;
        switch(mt){
            case 0:
                m.SetVentralMuscleInput(i, NMJ_Gain(i)*ventralHeadInputA);
                break;
            case 1:
                m.SetVentralMuscleInput(i, NMJ_Gain(i)*((ventralHeadInputA + ventralHeadInputP)/2));
                break;
            case 2:
                m.SetVentralMuscleInput(i, NMJ_Gain(i)*ventralHeadInputP);
                break;
        }
    }

    // Update Muscle activation
    m.EulerStep(StepSize);

    // Set input to Body
    //  First two segments receive special treatment because they are only affected by a single muscle
    b.SetDorsalSegmentActivation(1, m.DorsalMuscleOutput(1)/2);
    b.SetVentralSegmentActivation(1, m.VentralMuscleOutput(1)/2);
    b.SetDorsalSegmentActivation(2, m.DorsalMuscleOutput(1)/2);
    b.SetVentralSegmentActivation(2, m.VentralMuscleOutput(1)/2);

    //  All other segments receive force from two muscles
    for (int i = 3; i <= N_segments-2; i++)
    {
        mi = (int) ((i-1)/2);
        b.SetDorsalSegmentActivation(i, (m.DorsalMuscleOutput(mi) + m.DorsalMuscleOutput(mi+1))/2);
        b.SetVentralSegmentActivation(i, (m.VentralMuscleOutput(mi) + m.VentralMuscleOutput(mi+1))/2);
    }

    //  Last two segments receive special treatment because they are only affected by a single muscle
    b.SetDorsalSegmentActivation(N_segments-1, m.DorsalMuscleOutput(N_muscles)/2);
    b.SetVentralSegmentActivation(N_segments-1, m.VentralMuscleOutput(N_muscles)/2);
    b.SetDorsalSegmentActivation(N_segments, m.DorsalMuscleOutput(N_muscles)/2);
    b.SetVentralSegmentActivation(N_segments, m.VentralMuscleOutput(N_muscles)/2);

    // Time
    t += StepSize;
}

double Worm::CoMx()
{
    double temp = 0.0;
    for (int i = 1; i <= N_rods; i++) {
        temp += b.X(i);
    }
    return temp/N_rods;
}

double Worm::CoMy()
{
    double temp = 0.0;
    for (int i = 1; i <= N_rods; i++) {
        temp += b.Y(i);
    }
    return temp/N_rods;
}

void Worm::Curvature(TVector<double> &c)
{
    double dx1,dy1,dx2,dy2,a,a1,a2,seg;
    int k=1;

    for (int i = 3; i < N_segments-1; i+=2)
    {
        dx1 = b.X(i) - b.X(i-2);
        dy1 = b.Y(i) - b.Y(i-2);
        dx2 = b.X(i+2) - b.X(i);
        dy2 = b.Y(i+2) - b.Y(i);

        a1 = atan2(dy1,dx1);
        a2 = atan2(dy2,dx2);

        if (a1 > PI/2 and a2 < -PI/2)
            a = (a1 - 2*PI) - a2;
        else
            if (a1 < -PI/2 and a2 > PI/2)
                a = a1 - (a2 - 2*PI);
            else
                a = a1-a2;

        seg = sqrt(pow(b.X(i-2)-b.X(i+2),2) + pow(b.Y(i-2)-b.Y(i+2),2));
        c(k) = (2*sin(a)/seg)/1000;
        k++;
    }
}

double Worm::Orientation()
{
    return atan2(b.Y(Head)-b.Y(Tail),b.X(Head)-b.X(Tail));
}

// Dump the state to OFS if SKIPS steps have been performed

void Worm::DumpBodyState(ofstream &ofs, int skips)
{
    static int tt = skips;

    if (++tt >= skips) {
        tt = 0;

        ofs << t;
        // Body
        for (int i = 1; i <= N_rods; i++)
        {
            ofs <<  " " << b.X(i) << " " << b.Y(i) << " " << b.Phi(i);
        }
        ofs << "\n";
    }
}

void Worm::DumpActState(ofstream &ofs, int skips)
{
    static int tt = skips;

    if (++tt >= skips) {
        tt = 0;

        ofs << t;
        //ofs << "\nSR: ";
        // Stretch receptors
        ofs <<  " " << sr.HeadDorsalOutput() << " " << sr.HeadVentralOutput();
        for (int i = 1; i <= N_stretchrec; i++) {
            ofs <<  " " << sr.VCDorsalOutput(i) << " " << sr.VCVentralAOutput(i) << " " << sr.VCVentralPOutput(i);;
        }
        // Head Neurons
        //ofs << "\nH: ";
        for (int i = 1; i <= 4; i++) {
            ofs <<  " " << h.NeuronOutput(i);
        }
        // Ventral Cord Motor Neurons
        //ofs << "\nV: ";
        for (int i = 1; i <= N_units; i++) {
            for (int j = 1; j <= N_neuronsperunit; j++) {
                ofs <<  " " << n.NeuronOutput(nn(j,i));
            }
        }
        // Muscles
        //ofs << "\nM: ";
        for (int i = 1; i <= N_muscles; i++) {
            ofs <<  " " << m.DorsalMuscleOutput(i) << " " << m.VentralMuscleOutput(i);
        }
        ofs << "\n";
    }
}

void Worm::DumpVoltage(ofstream &ofs, int skips)
{
    static int tt = skips;

    if (++tt >= skips) {
        tt = 0;

        ofs << t;
        // Head Neurons
        for (int i = 1; i <= 4; i++) {
            ofs <<  " " << h.NeuronState(i);
        }
        // Ventral Cord Motor Neurons
        for (int i = 1; i <= N_units; i++) {
            for (int j = 1; j <= N_neuronsperunit; j++) {
                ofs <<  " " << n.NeuronState(nn(j,i));
            }
        }
        ofs << "\n";
    }
}

void Worm::DumpParams(ofstream &ofs)
{
    ofs << "Time-constants: \n DB: " << n.NeuronTimeConstant(DB) << "\n VBA/P: " << n.NeuronTimeConstant(VBA) << " / " << n.NeuronTimeConstant(VBP) << "\n DD: " << n.NeuronTimeConstant(DD) << "\n VDA/P: " << n.NeuronTimeConstant(VDA) << " / " << n.NeuronTimeConstant(VDP) << endl;
    ofs << "Biases: \n DB: " << n.NeuronBias(DB) << "\n VBA/P: " << n.NeuronBias(VBA) << " / " << n.NeuronBias(VBP)  <<  "\n DD: " << n.NeuronBias(DD) << "\n VDA/P: " << n.NeuronBias(VDA) <<  " / " << n.NeuronBias(VDP) << endl;
    ofs << "Self conns: \n DB: " << n.ChemicalSynapseWeight(DB, DB) << "\n VBA/P: " << n.ChemicalSynapseWeight(VBA, VBA) << " / " << n.ChemicalSynapseWeight(VBP, VBP) << "\n DD: " << n.ChemicalSynapseWeight(DD, DD) <<  "\n VDA/P: " << n.ChemicalSynapseWeight(VDA, VDA) <<  " / " << n.ChemicalSynapseWeight(VDP, VDP) << endl;
    ofs << "Chem Conns: \n DB->DD: " << n.ChemicalSynapseWeight(DB, DD) <<  "\n DB->VDA/VDP: " << n.ChemicalSynapseWeight(DB, VDA) << " / " << n.ChemicalSynapseWeight(DB, VDP) << "\n VBA/P->DD: " << n.ChemicalSynapseWeight(VBA, DD) << " / " << n.ChemicalSynapseWeight(VBP, DD) << "\n VBA/P->VDA/P: " << n.ChemicalSynapseWeight(VBA, VDA) << " / " << n.ChemicalSynapseWeight(VBP, VDP) << "\n VDA/P->VBA/P: " << n.ChemicalSynapseWeight(VDA, VBA) << " / " << n.ChemicalSynapseWeight(VDP, VBP) << "\n DD->VDA: " << n.ChemicalSynapseWeight(DD, VDA) <<endl;
    ofs << "Gap Juncs: \n DB-DB+1: " << n.ElectricalSynapseWeight(DB, DB+N_neuronsperunit) << "\n VBA-VBP / VBP-VBP+1: " << n.ElectricalSynapseWeight(VBA, VBP) << " / " << n.ElectricalSynapseWeight(VBP, VBA+N_neuronsperunit) << "\n VBP-DB+1: " << n.ElectricalSynapseWeight(VBP, DB+N_neuronsperunit) << "\n DD-VDA/P: " << n.ElectricalSynapseWeight(DD, VDA) << " / " << n.ElectricalSynapseWeight(DD, VDP) << "\n DD-DD+1: " << n.ElectricalSynapseWeight(DD, DD+N_neuronsperunit) << "\n VDA-VDP / VDP-VDP+1: " << n.ElectricalSynapseWeight(VDA, VDP) << " / " << n.ElectricalSynapseWeight(VDP, VDA+N_neuronsperunit) <<  endl;
    ofs << "SR Gain (VNC and Head): " << sr.SRvncgain << " " << sr.SRheadgain << endl;
    ofs << "NMJ weights: \n B: " << NMJ_DB << " " << NMJ_VBa << " " << NMJ_VBp << "\n D: " <<  NMJ_DD << " " << NMJ_VDa << " " << NMJ_VDp << endl;
    ofs << "Head: \nBiases: \n SMD(D/V): " << h.NeuronBias(SMDD) << " " << h.NeuronBias(SMDV) << "\n RMD(D/V): "<< h.NeuronBias(RMDD) << " "<< h.NeuronBias(RMDV) << endl;
    ofs << "Time-constants: \n SMD(D/V): " << h.NeuronTimeConstant(SMDD) << " " << h.NeuronTimeConstant(SMDV) << "\n RMD(D/V): " << h.NeuronTimeConstant(RMDD) << " " << h.NeuronTimeConstant(RMDV) << endl;
    ofs << "Self conns: \n SMD(D/V): " <<h.ChemicalSynapseWeight(SMDD,SMDD) << " " << h.ChemicalSynapseWeight(SMDV,SMDV) << "\n RMD(D/V): " << h.ChemicalSynapseWeight(RMDD,RMDD) << " "<< h.ChemicalSynapseWeight(RMDV,RMDV) << endl;
    ofs << "Chem conns: " << "\n SMDD->RMDD: " << h.ChemicalSynapseWeight(SMDD, RMDD) << "\n SMDD->SMDV: " << h.ChemicalSynapseWeight(SMDD, SMDV) << "\n SMDD->RMDV: " << h.ChemicalSynapseWeight(SMDD, RMDV) << "\n RMDD->SMDD: " << h.ChemicalSynapseWeight(RMDD, SMDD) << "\n RMDD->SMDV: " << h.ChemicalSynapseWeight(RMDD, SMDV) << "\n RMDD->RMDV: " << h.ChemicalSynapseWeight(RMDD, RMDV) << "\n SMDV->SMDD: " << h.ChemicalSynapseWeight(SMDV, SMDD) << "\n SMDV->RMDD: " << h.ChemicalSynapseWeight(SMDV, RMDD) << "\n SMDV->RMDV: " << h.ChemicalSynapseWeight(SMDV, RMDV) << "\n RMDV->SMDD: " << h.ChemicalSynapseWeight(RMDV, SMDD) << "\n RMDV->RMDD: " << h.ChemicalSynapseWeight(RMDV, RMDD) << "\n RMDV->SMDV: " << h.ChemicalSynapseWeight(RMDV, SMDV) << endl;
    ofs << "Gap Juncs: " << "\n SMD-RMD: " << h.ElectricalSynapseWeight(SMDD, RMDD)
                         << "\n RMD-RMD: " << h.ElectricalSynapseWeight(RMDD, RMDV) << endl;
    ofs << "NMJ weights: \n SMD(D/V): " << NMJ_SMDD << " " << NMJ_SMDV << "\n RMD(D/V): " <<  NMJ_RMDD << " " << NMJ_RMDV << endl;
    ofs << "NMJ Gain: " << NMJ_Gain_Map << endl;
}
