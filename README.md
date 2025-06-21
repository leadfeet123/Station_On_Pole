# Product Requirements Document: Compact In-Line Solid-State Substation (CSSS)

## 1. Introduction

### 1.1. Purpose
This document outlines the requirements for a comprehensive Simulink model of a **Compact In-Line, Non-Isolated Solid-State Transformer (SST)**. This model will serve as a virtual prototype to simulate, analyze, and validate the performance, control strategies, and operational limits of this novel substation concept before any physical prototyping.

### 1.2. Scope
The model will encompass the power electronics, the high-frequency autotransformer, and the associated control systems. It will simulate the conversion of **138kV** transmission-level AC voltage to standard distribution feeder voltages. The scope includes steady-state operation, dynamic response to load changes, and behavior under specific fault conditions, with a primary focus on the challenges introduced by the non-isolated architecture.

### 1.3. Objectives
* To create a functional and configurable model of the single-stage, non-isolated SST.
* To develop and test control algorithms for voltage regulation, power factor correction (PFC), and fault current limiting (FCL).
* To analyze system efficiency, power quality (harmonics), and component stresses under various operating conditions.
* To provide a platform for future research into protection schemes, grounding strategies, and system stability in a non-isolated configuration.

## 2. System Overview

### 2.1. Architecture
The model will represent a direct AC-AC SST based on the following "extreme decisions" for maximum compactness:

* **Single-Stage Conversion:** A direct AC-AC converter topology (e.g., Matrix or AC-AC Dual Active Bridge) will be used to minimize components and conversion losses.
* **Non-Isolated Autotransformer:** A high-frequency autotransformer with an amorphous magnetic core will provide voltage step-down *without* galvanic isolation.
* **In-Line Deployment Concept:** The model's parameters will be based on the assumption of a compact design mounted on a steel structure within a transmission right-of-way.

### 2.2. High-Level Diagram
```
      +-----------------+   +----------------------+   +-----------------+
      |   138kV Grid    |   |  Single-Stage AC-AC  |   |  Distribution   |
      |     Source      +--->|      Converter       +--->|      Load       |
      +-----------------+   | (with HF Autoxfmr)   |   +-----------------+
                          +----------------------+
                          ^
                          |
      +----------------------+
      |    Control System    |
      | (V-reg, PFC, FCL)    |
      +----------------------+
```

## 3. Functional Requirements

| ID      | Requirement                | Description                                                                                                                                              | Priority  |
| :------ | :------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------- | :-------- |
| `FR-01` | **Voltage Transformation** | The model shall accurately step down the input voltage from 138kV (RMS, line-to-line) to a configurable distribution voltage (e.g., 34.5kV or 11kV).     | Must-Have |
| `FR-02` | **Bidirectional Power Flow** | The model shall support power flow in both forward (transmission to distribution) and reverse (distribution to transmission) directions.                 | Must-Have |
| `FR-03` | **Output Voltage Regulation**| The control system shall maintain the output voltage at the distribution feeder within ±2% of the nominal setpoint under varying load conditions.        | Must-Have |
| `FR-04` | **Power Factor Correction** | The SST shall maintain a power factor of at least 0.99 (leading or lagging) at the 138kV input terminal.                                                 | Must-Have |
| `FR-05` | **Fault Current Limiting** | Upon detection of a downstream fault, the SST shall limit the fault current contribution to a configurable value (e.g., 1.5x rated current) within 2ms. | Must-Have |
| `FR-06` | **Harmonic Distortion** | The Total Harmonic Distortion (THD) of the output voltage waveform shall be maintained below 3% under nominal operating conditions.                      | Should-Have |
| `FR-07` | **Efficiency Analysis** | The model shall allow for the calculation of overall system efficiency by accounting for conduction, switching, and magnetic losses.                       | Should-Have |

## 4. Non-Functional Requirements

| ID       | Requirement                 | Description                                                                                                                                 | Priority  |
| :------- | :-------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------ | :-------- |
| `NFR-01` | **Parameterization** | All key component parameters (e.g., SiC device characteristics, core properties, control gains) shall be user-configurable from a single script. | Must-Have |
| `NFR-02` | **Modularity** | The model shall be structured modularly, with separate subsystems for the power stage, autotransformer, and control system.                   | Must-Have |
| `NFR-03` | **Simulation Performance** | The model should be optimized for a reasonable simulation speed, allowing for several seconds of operation to be simulated in minutes.          | Should-Have |
| `NFR-04` | **Documentation** | The Simulink model shall be well-documented with annotations explaining the function of key blocks and subsystems.                            | Could-Have |

## 5. Model Specifications

| Parameter                 | Value / Range                                                                                  |
| :------------------------ | :--------------------------------------------------------------------------------------------- |
| **Input Voltage** | 138 kV (AC, 3-phase, 60 Hz)                                                                    |
| **Output Voltage** | 11 kV - 34.5 kV (configurable)                                                                 |
| **Power Rating** | 10 MVA                                                                                         |
| **Converter Topology** | Direct AC-AC (e.g., Matrix or AC-AC DAB)                                                       |
| **Switching Frequency** | 10 kHz - 50 kHz (configurable)                                                                 |
| **Semiconductor Devices** | Silicon Carbide (SiC) MOSFETs (modeled with realistic on-resistance and switching characteristics) |
| **Autotransformer Core** | Amorphous Alloy (modeled with non-linear B-H curve for saturation and core loss calculations)    |

## 6. Assumptions and Dependencies

* The model will initially assume an ideal grid source, with the option to later introduce disturbances.
* The model will not include detailed thermal modeling but will calculate power losses that can be used as inputs for a separate thermal analysis.
* The model operates on the principle of **no galvanic isolation**. The safety and grounding implications will be analyzed through simulation results rather than being modeled as a separate safety system.
* The model will be developed using **MATLAB/Simulink R2024a or later** with the **Simscape Electrical™** toolbox.

## 7. Success Metrics

* **Primary Metric:** Successful, stable simulation of the CSSS under all scenarios defined in the Development Plan (steady-state, load changes, fault conditions).
* **Secondary Metrics:**
    * Achieving the target performance for voltage regulation (FR-03) and PFC (FR-04).
    * Demonstration of effective fault current limiting (FR-05).
    * Generation of a comprehensive report analyzing the trade-offs of the non-isolated design.
