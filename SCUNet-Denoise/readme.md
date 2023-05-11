Step 1: Download the SCUNet models
python main_download_pretrained_models.py --models "SCUNet" --model_dir "model_zoo"

Step 2: Add the input images into the folder testsets/input-img

Step 3: Image Denoising 
python main_test_scunet_real_application.py --model_name scunet_color_real_psnr --testset_name input-img