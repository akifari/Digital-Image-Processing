using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Drawing.Imaging;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;
using System.IO;


namespace dip
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }
        Bitmap SRC_IMG, DEST_IMG, Temp_IMG;


        private Bitmap ScaleImage(Bitmap image)
        {
            double ratioX = (double)pictureBox1.Width / image.Width;
            double ratioY = (double)pictureBox1.Height / image.Height;
            double ratio = Math.Min(ratioX, ratioY);
            int w = (int)(image.Width * ratio);
            int h = (int)(image.Height * ratio);
            Bitmap result = new Bitmap(image, new Size(w, h));
            return result;
        }
        private void button1_Click(object sender, EventArgs e)
        {
            openFileDialog1.ShowDialog();
            string File_Name = openFileDialog1.FileName.Trim();
            try
            {
                if (File_Name != "") pictureBox1.Load(File_Name);

                SRC_IMG = new Bitmap(File_Name);

                Bitmap display = ScaleImage(SRC_IMG);

                pictureBox1.Image = display;

                DEST_IMG = new Bitmap(SRC_IMG.Width, SRC_IMG.Height);
                Temp_IMG = new Bitmap(SRC_IMG.Width, SRC_IMG.Height);

            }

            catch
            {
                MessageBox.Show("Error Wrong File Type");
            }
        }





        private void button3_Click(object sender, EventArgs e)
        {
            int i, j;



            for (i = 0; i < SRC_IMG.Width; i++)
                for (j = 0; j < SRC_IMG.Height; j++)
                {
                    int RR = SRC_IMG.GetPixel(i, j).R;
                    int GG = SRC_IMG.GetPixel(i, j).G;
                    int BB = SRC_IMG.GetPixel(i, j).B;
                    byte GRAY = (byte)(.299 * RR + .587 * GG + .114 * BB);
                    DEST_IMG.SetPixel(i, j, Color.FromArgb(GRAY, GRAY, GRAY));
                }
            Bitmap display_2 = ScaleImage(DEST_IMG);
            pictureBox2.Image = display_2;
        }

        private void button2_Click(object sender, EventArgs e)
        {
            SRC_IMG = DEST_IMG;
            pictureBox1.Image = SRC_IMG;


        }



        private double[,] GaussianKernel(int lenght, double weight)
        {
            double[,] kernel = new double[lenght, lenght];
            double kernelSum = 0;
            int foff = (lenght - 1) / 2;
            double distance = 0;
            double constant = 1d / (2 * Math.PI * weight * weight);
            for (int y = -foff; y <= foff; y++)
            {
                for (int x = -foff; x <= foff; x++)
                {
                    distance = ((y * y) + (x * x)) / (2 * weight * weight);
                    kernel[y + foff, x + foff] = constant * Math.Exp(-distance);
                    kernelSum += kernel[y + foff, x + foff];
                }
            }
            for (int y = 0; y < lenght; y++)
            {
                for (int x = 0; x < lenght; x++)
                {
                    kernel[y, x] = kernel[y, x] * 1d / kernelSum;
                }
            }
            return kernel;
        }

        public static Bitmap Convolve(Bitmap srcImage, double[,] kernel)
        {
            int width = srcImage.Width;
            int height = srcImage.Height;
            BitmapData srcData = srcImage.LockBits(new Rectangle(0, 0, width, height),
                ImageLockMode.ReadOnly, PixelFormat.Format32bppArgb);
            int bytes = srcData.Stride * srcData.Height;
            byte[] buffer = new byte[bytes];
            byte[] result = new byte[bytes];
            Marshal.Copy(srcData.Scan0, buffer, 0, bytes);
            srcImage.UnlockBits(srcData);
            int colorChannels = 3;
            double[] rgb = new double[colorChannels];
            int foff = (kernel.GetLength(0) - 1) / 2;
            int kcenter = 0;
            int kpixel = 0;
            for (int y = foff; y < height - foff; y++)
            {
                for (int x = foff; x < width - foff; x++)
                {
                    for (int c = 0; c < colorChannels; c++)
                    {
                        rgb[c] = 0.0;
                    }
                    kcenter = y * srcData.Stride + x * 4;
                    for (int fy = -foff; fy <= foff; fy++)
                    {
                        for (int fx = -foff; fx <= foff; fx++)
                        {
                            kpixel = kcenter + fy * srcData.Stride + fx * 4;
                            for (int c = 0; c < colorChannels; c++)
                            {
                                rgb[c] += (double)(buffer[kpixel + c]) * kernel[fy + foff, fx + foff];
                            }
                        }
                    }
                    for (int c = 0; c < colorChannels; c++)
                    {
                        if (rgb[c] > 255)
                        {
                            rgb[c] = 255;
                        }
                        else if (rgb[c] < 0)
                        {
                            rgb[c] = 0;
                        }
                    }
                    for (int c = 0; c < colorChannels; c++)
                    {
                        result[kcenter + c] = (byte)rgb[c];
                    }
                    result[kcenter + 3] = 255;
                }
            }
            Bitmap resultImage = new Bitmap(width, height);
            BitmapData resultData = resultImage.LockBits(new Rectangle(0, 0, width, height),
                ImageLockMode.WriteOnly, PixelFormat.Format32bppArgb);
            Marshal.Copy(result, 0, resultData.Scan0, bytes);
            resultImage.UnlockBits(resultData);
            return resultImage;
        }

        private void openFileDialog1_FileOk(object sender, CancelEventArgs e)
        {

        }

        private void button5_Click(object sender, EventArgs e)
        {
            Bitmap transform = Convolve(SRC_IMG, GaussianKernel(5, 16d));
            Bitmap display = ScaleImage(transform);
            pictureBox2.Image = display;
        }

        public static Bitmap PadImage(Bitmap image, int padding)
        {
            int w = image.Width;
            int h = image.Height;

            BitmapData image_data = image.LockBits(
                new Rectangle(0, 0, w, h),
                ImageLockMode.ReadOnly,
                PixelFormat.Format24bppRgb);

            int bytes = image_data.Stride * image_data.Height;
            byte[] buffer = new byte[bytes];

            Marshal.Copy(image_data.Scan0, buffer, 0, bytes);
            image.UnlockBits(image_data);

            int wpad = w + 2 * padding;
            int hpad = h + 2 * padding;

            Bitmap padded = new Bitmap(wpad, hpad);
            BitmapData padded_data = padded.LockBits(
                new Rectangle(0, 0, wpad, hpad),
                ImageLockMode.WriteOnly,
                PixelFormat.Format24bppRgb);

            int padded_bytes = padded_data.Stride * padded_data.Height;
            byte[] result = new byte[padded_bytes];
            for (int i = 0; i < padded_bytes; i++)
            {
                result[i] = 1;
            }

            for (int x = 0; x < w; x++)
            {
                for (int y = 0; y < h; y++)
                {
                    int img_pixel = x * 3 + y * image_data.Stride;
                    int pad_pixel = x * 3 + y * padded_data.Stride;
                    for (int c = 0; c < 3; c++)
                    {
                        result[padded_data.Stride * padding + 3 * padding + pad_pixel + c] = buffer[img_pixel + c];
                    }
                }
            }

            Marshal.Copy(result, 0, padded_data.Scan0, padded_bytes);
            padded.UnlockBits(padded_data);

            return padded;
        }



        public static double median(byte[] array)
        {
            int len = array.Length;
            int mid_idx = len / 2;

            var sorted = array.OrderBy(n => n);

            double median;

            if ((len % 2) == 0)
            {
                median = (sorted.ElementAt(mid_idx) + sorted.ElementAt(mid_idx - 1)) / 2;
            }

            else
            {
                median = sorted.ElementAt(mid_idx);
            }

            return median;
        }

        public static Bitmap MedianFilter(Bitmap image)
        {
            int w = image.Width;
            int h = image.Height;


            BitmapData image_data = image.LockBits(
                new Rectangle(0, 0, w, h),
                ImageLockMode.ReadOnly,
                PixelFormat.Format24bppRgb);

            int bytes = image_data.Stride * image_data.Height;
            byte[] buffer = new byte[bytes];

            Marshal.Copy(image_data.Scan0, buffer, 0, bytes);
            image.UnlockBits(image_data);

            int r = 1;

            int wres = w - 2 * r;
            int hres = h - 2 * r;

            Bitmap result_image = new Bitmap(wres, hres);
            BitmapData result_data = result_image.LockBits(
                new Rectangle(0, 0, wres, hres),
                ImageLockMode.WriteOnly,
                PixelFormat.Format24bppRgb);
            int res_bytes = result_data.Stride * result_data.Height;
            byte[] result = new byte[res_bytes];

            for (int x = r; x < w - r; x++)
            {
                for (int y = r; y < h - r; y++)
                {

                    int pixel_location = x * 3 + y * image_data.Stride;
                    int res_pixel_loc = (x - r) * 3 + (y - r) * result_data.Stride;
                    double[] median = new double[3];
                    byte[][] neighborhood = new byte[3][];

                    for (int c = 0; c < 3; c++)
                    {
                        neighborhood[c] = new byte[(int)Math.Pow(2 * r + 1, 2)];
                        int added = 0;
                        for (int kx = -r; kx <= r; kx++)
                        {
                            for (int ky = -r; ky <= r; ky++)
                            {
                                int kernel_pixel = pixel_location + kx * 3 + ky * image_data.Stride;
                                neighborhood[c][added] = buffer[kernel_pixel + c];
                                added++;
                            }
                        }
                    }

                    for (int c = 0; c < 3; c++)
                    {
                        result[res_pixel_loc + c] = (byte)(Form1.median(neighborhood[c]));

                    }
                }
            }

            Marshal.Copy(result, 0, result_data.Scan0, res_bytes);
            result_image.UnlockBits(result_data);

            return result_image;
        }


        private void button4_Click(object sender, EventArgs e)
        {



            Bitmap processed = Form1.MedianFilter(SRC_IMG);
            Bitmap display = ScaleImage(processed);
            pictureBox2.Image = DEST_IMG;
            pictureBox2.Image = display;
            DEST_IMG = new Bitmap(SRC_IMG.Width, SRC_IMG.Height);


        }





        private void button6_Click(object sender, EventArgs e)
        {

            Bitmap processed = ArithmeticMean(SRC_IMG);
            Bitmap display = ScaleImage(processed);
            pictureBox2.Image = DEST_IMG;
            pictureBox2.Image = display;
            DEST_IMG = new Bitmap(SRC_IMG.Width, SRC_IMG.Height);
        }


        public static Bitmap ArithmeticMean(Bitmap image)
        {
            int w = image.Width;
            int h = image.Height;

            BitmapData image_data = image.LockBits(
                new Rectangle(0, 0, w, h),
                ImageLockMode.ReadOnly,
                PixelFormat.Format24bppRgb);

            int bytes = image_data.Stride * image_data.Height;
            byte[] buffer = new byte[bytes];

            Marshal.Copy(image_data.Scan0, buffer, 0, bytes);
            image.UnlockBits(image_data);

            int r = 1;

            int wres = w - 2 * r;
            int hres = h - 2 * r;

            Bitmap result_image = new Bitmap(wres, hres);
            BitmapData result_data = result_image.LockBits(
                new Rectangle(0, 0, wres, hres),
                ImageLockMode.WriteOnly,
                PixelFormat.Format24bppRgb);
            int res_bytes = result_data.Stride * result_data.Height;
            byte[] result = new byte[res_bytes];


            for (int x = r; x < w - r; x++)
            {
                for (int y = r; y < h - r; y++)
                {

                    int pixel_location = x * 3 + y * image_data.Stride;
                    int res_pixel_loc = (x - r) * 3 + (y - r) * result_data.Stride;
                    double[] mean = new double[3];

                    for (int kx = -r; kx <= r; kx++)
                    {
                        for (int ky = -r; ky <= r; ky++)
                        {

                            int kernel_pixel = pixel_location + kx * 3 + ky * image_data.Stride;

                            for (int c = 0; c < 3; c++)
                            {
                                mean[c] += buffer[kernel_pixel + c] / Math.Pow(2 * r + 1, 2);
                            }
                        }
                    }

                    for (int c = 0; c < 3; c++)
                    {
                        result[res_pixel_loc + c] = (byte)mean[c];
                    }
                }
            }

            Marshal.Copy(result, 0, result_data.Scan0, res_bytes);
            result_image.UnlockBits(result_data);

            return result_image;
        }

        private void button7_Click(object sender, EventArgs e)
        {

            Bitmap processed = ImageSharpen(SRC_IMG);
            Bitmap display = ScaleImage(processed);
            pictureBox2.Image = DEST_IMG;
            pictureBox2.Image = display;
            DEST_IMG = new Bitmap(SRC_IMG.Width, SRC_IMG.Height);

        }

        private void button8_Click(object sender, EventArgs e)
        {




            Bitmap processed = EdgeDetect(SRC_IMG);
            Bitmap display = ScaleImage(processed);
            pictureBox2.Image = DEST_IMG;
            pictureBox2.Image = display;
            
            DEST_IMG = new Bitmap(SRC_IMG.Width, SRC_IMG.Height);

            
        }

        private void pictureBox1_Click(object sender, EventArgs e)
        {

        }



        public static Bitmap ImageSharpen(Bitmap image)
        {
            int w = image.Width;
            int h = image.Height;

            BitmapData image_data = image.LockBits(
                new Rectangle(0, 0, w, h),
                ImageLockMode.ReadOnly,
                PixelFormat.Format24bppRgb);

            int bytes = image_data.Stride * image_data.Height;
            byte[] buffer = new byte[bytes];
            byte[] result = new byte[bytes];

            Marshal.Copy(image_data.Scan0, buffer, 0, bytes);
            image.UnlockBits(image_data);

            for (int i = 2; i < w - 2; i++)
            {
                for (int j = 2; j < h - 2; j++)
                {
                    int p = i * 3 + j * image_data.Stride;
                    for (int k = 0; k < 3; k++)
                    {
                        double val = 0d;
                        for (int xkernel = -1; xkernel < 2; xkernel++)
                        {
                            for (int ykernel = -1; ykernel < 2; ykernel++)
                            {
                                int kernel_p = k + p + xkernel * 3 + ykernel * image_data.Stride;
                                val += buffer[kernel_p] * Kernels.Laplacian[xkernel + 1, ykernel + 1];
                            }
                        }
                        val = val > 0 ? val : 0;
                        result[p + k] = (byte)((val + buffer[p + k]) > 255 ? 255 : (val + buffer[p + k]));
                    }
                }
            }

            Bitmap res_img = new Bitmap(w, h);
            BitmapData res_data = res_img.LockBits(
                new Rectangle(0, 0, w, h),
                ImageLockMode.WriteOnly,
                PixelFormat.Format24bppRgb);
            Marshal.Copy(result, 0, res_data.Scan0, bytes);
            res_img.UnlockBits(res_data);
            return res_img;
        }


        public static class Kernels
        {
            public static double[,] Laplacian
            {
                get
                {
                    return new double[,]
                    {
                    { 0,-1, 0 },
                    {-1, 4,-1 },
                    { 0,-1, 0 }
                    };
                }
            }
        }

        public static byte[] Convolute(byte[] buffer, BitmapData image_data, int[,] filter)
        {
            byte[] result = new byte[buffer.Length];
            int ho = (filter.GetLength(0) - 1) / 2;
            int vo = (filter.GetLength(1) - 1) / 2;

            for (int x = ho; x < image_data.Width - ho; x++)
            {
                for (int y = vo; y < image_data.Height - vo; y++)
                {
                    int position = x * 3 + y * image_data.Stride;
                    int sum = 0;
                    for (int i = -ho; i <= ho; i++)
                    {
                        for (int j = -vo; j <= vo; j++)
                        {
                            int filter_position = position + i * 3 + j * image_data.Stride;
                            sum += (buffer[filter_position] * filter[i + ho, j + vo]);
                        }
                    }

                    for (int c = 0; c < 3; c++)
                    {
                        if (sum > 255)
                        {
                            sum = 255;
                        }

                        else if (sum < 0)
                        {
                            sum = Math.Abs(sum);
                        }

                        result[position + c] = (byte)(sum);
                    }
                }
            }

            return result;
        }


        


        private void button9_Click(object sender, EventArgs e)
        {
           
            if (pictureBox1.Image != null)
            {

                Bitmap processed = new Bitmap(pictureBox1.Image);
                processed.RotateFlip(RotateFlipType.Rotate180FlipX);
                Bitmap display = ScaleImage(processed);
                pictureBox2.Image = DEST_IMG;
                pictureBox2.Image = display;
                
            


            }
        }

        private void button10_Click(object sender, EventArgs e)
        {

            if (pictureBox1.Image != null)
            {

                Bitmap processed = new Bitmap(pictureBox1.Image);
                processed.RotateFlip(RotateFlipType.Rotate90FlipXY);
                Bitmap display = ScaleImage(processed);
                pictureBox2.Image = DEST_IMG;
                pictureBox2.Image = display;




            }
        }

        public static Bitmap EdgeDetect(Bitmap image)
        {
            int w = image.Width;
            int h = image.Height;

            BitmapData image_data = image.LockBits(
                new Rectangle(0, 0, w, h),
                ImageLockMode.ReadOnly,
                PixelFormat.Format24bppRgb);

            int bytes = image_data.Stride * image_data.Height;
            byte[] buffer = new byte[bytes];
            byte[] temp = new byte[bytes];
            byte[] result = new byte[bytes];

            Marshal.Copy(image_data.Scan0, buffer, 0, bytes);
            image.UnlockBits(image_data);

            for (int i = 0; i < bytes; i++)
            {
                result[i] = (byte)((result[i] + temp[i]) > 255 ? 255 : (result[i] + temp[i]));
            }

            Bitmap res_img = new Bitmap(w, h);
            BitmapData res_data = res_img.LockBits(
                new Rectangle(0, 0, w, h),
                ImageLockMode.WriteOnly,
                PixelFormat.Format24bppRgb);
            Marshal.Copy(result, 0, res_data.Scan0, bytes);
            res_img.UnlockBits(res_data);

            return res_img;
        }

        public static class Filters
        {
            public static int[,] PrewittVertical
            {
                get
                {
                    return new int[,]
                    {
                    { -1, -1, -1 },
                    { 0, 0, 0 },
                    { 1, 1, 1 }
                    };
                }
            }

            public static int[,] PrewittHorizontal
            {
                get
                {
                    return new int[,]
                    {
                    { -1, 0, 1 },
                    { -1, 0, 1 },
                    { -1, 0, 1 }
                    };
                }
            }










        }














    }
}

