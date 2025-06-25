# Content-Aware Image Resizing (Seam Carving with Adaptive Energy Maps)

This project presents an advanced approach to image resizing using seam carving, enhanced with a novel adaptive energy map that combines multiple cues such as depth, saliency, gradient, and edges. The goal is to minimize visual distortion and preserve important content when reducing image size.

---

## Introduction

In modern applications like web design, social media, and mobile apps, efficient image resizing is essential. Traditional methods like bicubic interpolation often distort or remove important content. Seam carving is a content-aware resizing technique that overcomes this by identifying and removing the least noticeable paths (seams) in an image.

Our project builds upon this by introducing an adaptive energy map, using a linear combination of several energy cues with carefully tuned weights.

---

## Method Overview

### Energy Map Composition

We generate five different energy maps:

- **Depth map**
- **Saliency map**
- **Gradient magnitude**
- **Sobel edge detection**
- **Canny edge detection**

Each energy map is assigned a custom weight based on its effectiveness:

| Energy Type     | Weight |
|------------------|--------|
| Depth            | 3      |
| Gradient/Sobel/Canny (Edges) | 2 |
| Saliency         | 1      |

These are linearly combined to form the final **adaptive energy map**:

E_total = 3 * E_depth + 2 * (E_gradient + E_sobel + E_canny) + 1 * E_saliency


### Adaptive Resizing Strategy

- If the image has **high depth complexity**, we apply seam carving only to a portion of the image, and perform the rest of the resizing using **bicubic interpolation** to avoid losing critical structures.
- Otherwise, the full resize is performed using seam carving guided by the combined energy map.

### Reference Work

We were inspired by [Improved seam carving combining with 3D saliency for image retargeting](https://www.sciencedirect.com/science/article/abs/pii/S0925231214013769), where the authors combine depth and saliency adaptively based on normalized variance of the depth map:

E_s = (1 - α) * E_saliency + α * E_depth


---

## Results

We tested our method on several images with different structures and compared:

1. Seam carving using only **depth**
2. Seam carving using only **saliency**
3. Equal-weight linear combination of depth and saliency
4. **Our adaptive method** with full combination of energy maps

Our method showed better preservation of key image content, especially in complex scenes with depth variations.

![image](https://github.com/user-attachments/assets/3f5fb0d3-62e7-4181-87a1-978e4975db6c)
![image](https://github.com/user-attachments/assets/92512cde-5089-4f14-ab3a-8aa9411fa0f6)
![image](https://github.com/user-attachments/assets/07bb3c3c-3a89-4b17-a70a-b61c9ee83e1a)
![image](https://github.com/user-attachments/assets/3e04dc9e-5e04-4ee5-b0ff-aff41ce9623b)


