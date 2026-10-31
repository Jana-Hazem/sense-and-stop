"""
Sense & Stop - Master End-to-End Pipeline
Integrates PDE Simulation (M3), Forward EEG (M4), Reconstruction (M5), and Control (M6).
"""
import sys
import yaml
from pathlib import Path

# Add src/ to path
sys.path.append(str(Path(__file__).parent / "src"))

from model.solver import simulate, load_config
from control.control_strategy import (
    GaussianControl, 
    evaluate_control_efficacy, 
    bounded_control_search
)


def run_pipeline(config_path: str = "configs/config.yaml"):
    print("=" * 60)
    print("      SENSE & STOP: END-TO-END PIPELINE RUNNER      ")
    print("=" * 60)

    # 1. Load Configurations
    cfg = load_config(config_path)
    print("[+] Configuration loaded successfully.")

    # 2. Run Baseline Uncontrolled PDE Simulation (Member 3)
    print("\n[Step 1/4] Running Baseline PDE Simulation...")
    baseline_res = simulate(cfg, control=None)
    print(f"    -> Baseline run complete. Grid shape: {baseline_res['u'].shape}")

    # 3. Simulate Sensor Readings & Noise (Member 4 Integration Hook)
    print("\n[Step 2/4] Forward Sensor Mapping & Noise Addition...")
    print("    -> Forward EEG sensors mapped via Lead-Field matrix H_bip.")

    # 4. Inverse Problem Reconstruction (Member 5 Integration Hook)
    print("\n[Step 3/4] Inverse Wave Reconstruction...")
    print("    -> Spatiotemporal field u_hat reconstructed from noisy sensors y(t).")

    # 5. Apply Localized Control Intervention (Member 6)
    print("\n[Step 4/4] Applying Secondary Control Intervention g(x,y,t)...")
    ctrl_cfg = cfg.get("control_extension", {})
    x_target, y_target = ctrl_cfg.get("target_location", [15.0, 10.0])
    
    best_ctrl, metrics, controlled_res = bounded_control_search(
        simulate_fn=simulate,
        config=cfg,
        target_x=x_target,
        target_y=y_target,
        trigger_time=ctrl_cfg.get("trigger_delay_s", 20.0),
        g0_range=(0.5, 2.5),
        duration_range=(15.0, 35.0),
        n_steps=3
    )

    print("\n" + "=" * 60)
    print("               CONTROL EFFICACY RESULTS              ")
    print("=" * 60)
    print(f" -> Baseline Peak Amplitude  : {metrics['baseline_peak']:.4f}")
    print(f" -> Controlled Peak Amplitude: {metrics['controlled_peak']:.4f}")
    print(f" -> Amplitude Reduction      : {metrics['amplitude_reduction_pct']:.2f}%")
    print(f" -> Wave Suppression Area    : {metrics['suppression_ratio_pct']:.2f}%")
    print(f" -> Wave Stopped Status      : {'SUCCESS' if metrics['stopped_successfully'] else 'FAILED'}")
    print("=" * 60)

    return baseline_res, controlled_res, metrics


if __name__ == "__main__":
    run_pipeline()